#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<time.h>
#include <math.h>
#include<omp.h>

#include "FreeImage.h"
#include "mpi.h"

#define H 200
#define W 400
#define ITER 100//number of iterations
#define P_ST 0.01 //prob. of fire starting in a cell at init
#define P_BD 0.6 //prob. of burning cell burning down
#define P_H 0.58 //constant spreading prob.
#define IMAGE_PATH "./temp/" //image sequence is stored here

#define X_EDGE 1 //width of the horizontal "edge" of neighbouring cells to exhange
#define Y_EDGE 0 //width of the vertical edge (has to be 0 or same as X_edge)

//int grid[H][W];
//int next_grid[H][W];
int t;
unsigned int proc_seed;
int loc_h;
int loc_w;
char symbols[5] = {'E',' ','.','&','*'};
unsigned char *image;
int colors[5][3] = {{0,0,0},
                {80, 80, 80},
                {10, 40, 15},
                {255,215,80},
                {30, 30, 30 }};

/*  
    Cell states:
    1 - non-flammable cell
    2 - flammable but intact
    3 - currently burning
    4 - had burned down
*/

double rnd(){
    return (double) rand_r(&proc_seed) / (double) RAND_MAX;
}

// Store row to out_buf for send
void get_row(char** grid, char* out_buf, int y) {
    for (int x = 0; x < loc_w; x++){
        out_buf[x] = grid[y][x];
    }
}

// Store col to out_buf for send
void get_col(char** grid, char* out_buf, int x){
    for (int y = 0; y < loc_h; y++){
        out_buf[y] = grid[y][x];
    }
}

// Store row in in_buf to grid
void set_row(char** grid, char* in_buf, int y) {
    for (int x = 0; x < loc_w; x++){
        grid[y][x] = in_buf[x];
    }
}

// Store col in in_buf to grid
void set_col(char** grid, char* in_buf, int x) {
    for (int y = 0; y < loc_h; y++){
        grid[y][x] = in_buf[y];
    }
} 

// Initialize cells (to 1 or 2)
void init(char** grid, double density){
    for(int y = 0; y < loc_h; y++)
        for(int x = 0; x < loc_w; x++){
            grid[y][x] = rnd() < density ? 2 : 1;
        }
}

// Start the fire in some cells: n - number of sparks
void spark(char** grid, int n){
    int y0, x0;
    for(int i = 0; i < n; i++){
        y0 = rand_r(&proc_seed) % loc_h;
        x0 = rand_r(&proc_seed) % loc_w;
        grid[y0][x0] = 3;
    }
}

// Get state of a cell at t+1
int next_state(char** grid, int y, int x){
    if (grid[y][x] == 1)
        return 1;
    else if (grid[y][x] == 4)
        return 4;
    else if (grid[y][x] == 3 && rnd() < P_BD)
        return 4;
    else if (grid[y][x] == 2){
        for(int i = -1; i <= 1; i++)
            for(int j = -1; j <=1; j++){           
                // fire spreads from each neighbor with some probability
                double p_spread = P_H; //change to function
                if (y+i < 0-X_EDGE || y+i > loc_h+X_EDGE || x+j < 0-Y_EDGE || x+j > loc_w+Y_EDGE)
                    continue;
                if (grid[y+i][x+j] == 3 && rnd() < p_spread)
                    return 3;
            }

    }

    return grid[y][x];
}

// Print the grid in current state
void print_grid(char** grid, int* target_cells){
    for(int y = 0; y < loc_h; y++){
        for(int x = 0; x < loc_w; x++){
            printf("%c", symbols[grid[y][x]]);
        }
        printf("\n");
    }

    if (target_cells != NULL && 0 <= target_cells[0]) {
        printf("\nTOP\n");
        for (int i = 0; i < loc_w; i++){
            printf("%c", symbols[grid[-1][i]]);
        }
    }
    if (target_cells != NULL && 0 <= target_cells[1]) {
        printf("\nBOT\n");
        for (int i = 0; i < loc_w; i++){
            printf("%c", symbols[grid[loc_h][i]]);
        }
    }
         
}

// Parallel state update for whole grid with no send
void update_nosend(char** grid, char** next_grid){ 
    #pragma omp parallel
    {
        #pragma omp for collapse(2) schedule(static)
        for(int y = 1; y < loc_h-1; y++)
            for(int x = 1; x < loc_w-1; x++)
                next_grid[y][x] = next_state(grid,y,x);

        #pragma omp for collapse(2) schedule(static)
        for(int y = 1; y < loc_h-1; y++)
            for(int x = 1; x < loc_w-1; x++)
                grid[y][x] = next_grid[y][x];
    }
}

// Parallel state update with send functions
void update_send(char** grid, char** next_grid, int* target_cells, MPI_Comm grid_comm, int my_id) {

    MPI_Request send_up, send_down, recv_up, recv_down;
    MPI_Status status;

    // Prepare to send and receive asynchronous
    char* out_buf_top = (char*)malloc(loc_w * sizeof(char));
    char* out_buf_bot = (char*)malloc(loc_w * sizeof(char));
    char* in_buf_top = (char*)malloc(loc_w * sizeof(char));
    char* in_buf_bot = (char*)malloc(loc_w * sizeof(char));

    // Check for top neigbour and send upwards
    if (0 <= target_cells[0]) {
        get_row(grid, out_buf_top, 0);
        MPI_Isend(out_buf_top, loc_w, MPI_CHAR, target_cells[0], my_id, grid_comm, &send_up);
        MPI_Irecv(in_buf_top, loc_w, MPI_CHAR, target_cells[0], target_cells[0], grid_comm, &recv_up);
    }
    // Check for bot neigbour
    if (0 <= target_cells[1]) {
        get_row(grid, out_buf_bot, loc_h-1);
        MPI_Isend(out_buf_bot, loc_w, MPI_CHAR, target_cells[1], my_id, grid_comm, &send_down);
        MPI_Irecv(in_buf_bot, loc_w, MPI_CHAR, target_cells[1], target_cells[1], grid_comm, &recv_down);
    }

    // Calculate all except borders
    #pragma omp parallel
    {
        #pragma omp for collapse(2) schedule(guided)
        for(int y = 1; y < loc_h-1; y++)
            for(int x = 0; x < loc_w; x++)
                next_grid[y][x] = next_state(grid,y,x);
    }

    // Wait for send to top
    if (0 <= target_cells[0]) {
        MPI_Wait(&send_up, MPI_STATUS_IGNORE);
    }
    // Wait for recv from bot and send to bot
    if (0 <= target_cells[1]) {
        MPI_Wait(&recv_down, MPI_STATUS_IGNORE);
        set_row(grid, in_buf_bot, loc_h);
        MPI_Wait(&send_down, MPI_STATUS_IGNORE);
    }
    // Wait for recv from top
    if (0 <= target_cells[0]) {
        MPI_Wait(&recv_up, MPI_STATUS_IGNORE);
        set_row(grid, in_buf_top, -1);
    }

    // Calculate all borders and update grid
    #pragma omp parallel
    {
        #pragma omp for schedule(guided)
        for(int x = 0; x < loc_w; x++) {
            next_grid[0][x] = next_state(grid,0,x);
            next_grid[loc_h-1][x] = next_state(grid,loc_h-1,x); 
        }

        #pragma omp for collapse(2) schedule(guided)
        for(int y = 0; y < loc_h; y++)
            for(int x = 0; x < loc_w; x++)
                grid[y][x] = next_grid[y][x];
    }

    MPI_Barrier(grid_comm);
    free(out_buf_top);
    free(out_buf_bot);
    free(in_buf_top);
    free(in_buf_bot);
}

// Returns grid that can be accessed like [i][j], and has "out-of-bound" edges allocated
char** alloc_grid(int h, int w, int x_edge, int y_edge){
    char* alloc = (char *) malloc((h+2*x_edge) * (w+2*y_edge) * sizeof(char));
    char** grid = (char **) malloc((h+2*x_edge) * sizeof(char *));

    for(int i = 0; i < (h+2*x_edge); i++)
            grid[i] = &alloc[i*(w+2*y_edge)] + y_edge;
    grid += x_edge;

    return grid;
}

int main(int argc, char* argv[]){

    int num_procs, my_id;
    MPI_Status status;

    // Config variables for Cartesian
    int dim_sizes[2];
    int wrap_around[2];
    int dims;
    MPI_Comm grid_comm;

    int my_grid_id;
    int my_grid_coords[2];
    int my_neigbours[2];

    // ================== MPI ==================
    double t_start = omp_get_wtime();
    // Configure MPI parallelization based on input arguments
    MPI_Init(&argc, &argv);    
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    dims = 1;
    // If num of procs nod divisble by num of dimensions
    if (num_procs % dims != 0) {
        printf("Exit:%d %d %d\n", num_procs, dims, num_procs % dims);
        MPI_Abort(MPI_COMM_WORLD, 1);
        exit(1);
    }
    
    dim_sizes[0] = num_procs/dims; dim_sizes[1] = dims;
    wrap_around[0] = 0; wrap_around[1] = 0; // Set non circular cartesian grid

    MPI_Cart_create(MPI_COMM_WORLD, dims, dim_sizes, wrap_around, 0, &grid_comm);

    //Get coordinates
    MPI_Comm_rank(grid_comm, &my_grid_id);
    MPI_Cart_get(grid_comm, dims, dim_sizes, wrap_around, my_grid_coords);

    //Get neighbouring processes (for row split)
    MPI_Cart_shift(grid_comm, 0, 1, &my_neigbours[0], &my_neigbours[1]);

    // Each process generates its own grid
    int h_size = (int)round((double) H / (double) num_procs);
    if (my_id == num_procs-1) {
        //last process gets less
        h_size = H - h_size * (num_procs-1);
    }

    if (my_grid_id == 0) {
        printf("Res:%dx%d | grid: %dx%d\n", H, W, dim_sizes[0], dim_sizes[1]);
    }

    char** grid = alloc_grid(h_size, W, X_EDGE, Y_EDGE);
    char** next_grid = alloc_grid(h_size, W, X_EDGE, Y_EDGE);
    //printf("hsize: %d\n", h_size);
    
    // Setup local sizes for process
    loc_h = h_size;
    loc_w = W;    
    proc_seed = my_id;
    printf("id:%d [%d,%d] | neighb: [%d,%d] | sizes: [%d,%d]\n", my_grid_id, my_grid_coords[0], my_grid_coords[1], my_neigbours[0], my_neigbours[1], loc_h, loc_w);
    
    
    init(grid, 0.7);

    // Spark fire in three processes
    if (my_grid_id == 0)
        spark(grid, 5);

    /**
     * Simulation
     * 1. Create asynchronous send and receive
     * 2. Calculate all cells except borders
     * 3. Wait for send and recv
     * 4. Calculate borders and update grid
    */
    for (t = 0; t < ITER; t++) {
        update_send(grid,next_grid,my_neigbours,grid_comm,my_grid_id);
    }

    MPI_Finalize();

    double t_end = omp_get_wtime();
    if (my_grid_id == 0)
        printf("Elapsed time = %f\n", t_end-t_start);
}


// Compile with: mpicc -lm -fopenmp parallel.c -o out
// Run with: sbatch --reservation=fri ./job.sh