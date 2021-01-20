#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<time.h>
#include <math.h>
#include<omp.h>

#include "FreeImage.h"
#include "mpi.h"

#define H 8000
#define W 8000
#define ITER 1000//number of iterations
#define P_ST 0.01 //prob. of fire starting in a cell at init
#define P_BD 0.6 //prob. of burning cell burning down
#define P_H 0.58 //constant spreading prob.
#define IMAGE_PATH "./temp/" //image sequence is stored here
#define EDGE_LEN 1 //width of the edge (for row split horizontal and for block split both
#define BLOCKS 1 //0 - aggregate by rows and 1 - aggregate by blocks

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

// Find divisor closest to square root
int get_divisor(int x) {
    for (int j = (int)floor(sqrt(x)); j < x; j++){
        if (x % j == 0) 
            return j;
    }
    return 1;
}

// Get dimensions for Cartesian topology
int set_dims(int num_procs, int* dims) {
    if (BLOCKS) {
        int width = get_divisor(num_procs);
        if (width != 1) {
            dims[0] = num_procs/width;dims[1]=width;
            return 2;
        }
    }

    // Fallback use rows
    dims[0] = num_procs;dims[1]=1;
    return 1;
}

double rnd(){
    return (double) rand_r(&proc_seed) / (double) RAND_MAX;
}

// Store rows to out_buf for send from x to x+EDGE_LEN
void get_rows(char** grid, char* out_buf, int y) {
    int count = 0;
    while (count < EDGE_LEN) {
        for (int x = 0; x < loc_w; x++){
            out_buf[count*loc_w + x] = grid[y+count][x];
        }
        count++;
    }
}

// Store cols to out_buf for send from y to y+EDGE_LEN
void get_cols(char** grid, char* out_buf, int x){
    int count = 0;
    while (count < EDGE_LEN) {
        for (int y = 0; y < loc_h; y++){
            out_buf[count*loc_h + y] = grid[y][x + count];
        }
        count++;
    }
}

// Store rows in in_buf to grid from x to x+EDGE_LEN
void set_rows(char** grid, char* in_buf, int y) {
    int count = 0;
    while (count < EDGE_LEN) {
        for (int x = 0; x < loc_w; x++){
            grid[y+count][x] = in_buf[x + count*loc_w];
        }
        count++;
    }
}

// Store cols in in_buf to grid y to y+EDGE_LEN
void set_cols(char** grid, char* in_buf, int x) {
    int count = 0;
    while (count < EDGE_LEN) {
        for (int y = 0; y < loc_h; y++){
            grid[y][x+count] = in_buf[y + count*loc_h];
        }
        count++;
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
                if (y+i < 0-EDGE_LEN || y+i > loc_h+EDGE_LEN || x+j < 0-EDGE_LEN || x+j > loc_w+EDGE_LEN)
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
void update_nosend(char** grid, char** next_grid, int exchange, int blocks){ 
    int eH = - EDGE_LEN + 1 + exchange;
    int eW = blocks ? -EDGE_LEN + 1 + exchange : 0; 
    #pragma omp parallel
    {
        #pragma omp for collapse(2) schedule(guided)
        for(int y = eH; y < loc_h-eH; y++)
            for(int x = eW; x < loc_w-eW; x++)
                next_grid[y][x] = next_state(grid,y,x);

        #pragma omp for collapse(2) schedule(guided)
        for(int y = eH; y < loc_h-eH; y++)
            for(int x = eW; x < loc_w-eW; x++)
                grid[y][x] = next_grid[y][x];
    }
}

// Parallel state update with send functions
void update_send(char** grid, char** next_grid, int* target_cells, MPI_Comm grid_comm, int my_id) {

    MPI_Request send_up, send_down, recv_up, recv_down;
    MPI_Status status;

    // Prepare to send and receive asynchronous
    char* out_buf_top = (char*)malloc(EDGE_LEN * loc_w * sizeof(char));
    char* out_buf_bot = (char*)malloc(EDGE_LEN * loc_w * sizeof(char));
    char* in_buf_top = (char*)malloc(EDGE_LEN * loc_w * sizeof(char));
    char* in_buf_bot = (char*)malloc(EDGE_LEN * loc_w * sizeof(char));

    // Check for top neigbour and send upwards
    if (0 <= target_cells[0]) {
        get_rows(grid, out_buf_top, 0);
        MPI_Isend(out_buf_top, EDGE_LEN * loc_w, MPI_CHAR, target_cells[0], my_id, grid_comm, &send_up);
        MPI_Irecv(in_buf_top, EDGE_LEN * loc_w, MPI_CHAR, target_cells[0], target_cells[0], grid_comm, &recv_up);
    }

    // Check for bot neigbour
    if (0 <= target_cells[1]) {
        get_rows(grid, out_buf_bot, loc_h-EDGE_LEN);
        MPI_Isend(out_buf_bot, EDGE_LEN * loc_w, MPI_CHAR, target_cells[1], my_id, grid_comm, &send_down);
        MPI_Irecv(in_buf_bot, EDGE_LEN * loc_w, MPI_CHAR, target_cells[1], target_cells[1], grid_comm, &recv_down);
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
        set_rows(grid, in_buf_bot, loc_h);
        MPI_Wait(&send_down, MPI_STATUS_IGNORE);
    }
    // Wait for recv from top
    if (0 <= target_cells[0]) {
        MPI_Wait(&recv_up, MPI_STATUS_IGNORE);
        set_rows(grid, in_buf_top, -EDGE_LEN);
    }

    // Calculate all borders and update grid
    #pragma omp parallel
    {
        #pragma omp for collapse(2) schedule(guided)
        for(int x = 0; x < loc_w; x++) {
            for(int count = 0; count < EDGE_LEN; count++) {
                next_grid[0-count][x] = next_state(grid,0-count,x);
                next_grid[loc_h-1+count][x] = next_state(grid,loc_h-1+count,x);
            }
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

// Parallel state update with send functions for blocks
void update_send_blocks(char** grid, char** next_grid, int* target_cells, MPI_Comm grid_comm, int my_id){

    MPI_Request top_send,top_recv,bot_send,bot_recv,left_send,left_recv,right_send,right_recv;
    MPI_Status status;

    char* out_buf_top = (char*)malloc(EDGE_LEN * loc_w * sizeof(char));
    char* out_buf_bot = (char*)malloc(EDGE_LEN * loc_w * sizeof(char));
    char* out_buf_left = (char*)malloc(EDGE_LEN * loc_h * sizeof(char));
    char* out_buf_right = (char*)malloc(EDGE_LEN * loc_h * sizeof(char));

    char* in_buf_top = (char*)malloc(EDGE_LEN * loc_w * sizeof(char));
    char* in_buf_bot = (char*)malloc(EDGE_LEN * loc_w * sizeof(char));
    char* in_buf_left = (char*)malloc(EDGE_LEN * loc_h * sizeof(char));
    char* in_buf_right = (char*)malloc(EDGE_LEN * loc_h * sizeof(char));
    
    // Top
    if (0 <= target_cells[0]) {
        get_rows(grid, out_buf_top, 0);
        MPI_Isend(out_buf_top, EDGE_LEN * loc_w, MPI_CHAR, target_cells[0], my_id, grid_comm, &top_send);
        MPI_Irecv(in_buf_top, EDGE_LEN * loc_w, MPI_CHAR, target_cells[0], target_cells[0], grid_comm, &top_recv);
    }
    
    // Bot
    if (0 <= target_cells[1]) {
        get_rows(grid, out_buf_bot, loc_h-EDGE_LEN);
        MPI_Isend(out_buf_bot, EDGE_LEN * loc_w, MPI_CHAR, target_cells[1], my_id, grid_comm, &bot_send);
        MPI_Irecv(in_buf_bot, EDGE_LEN * loc_w, MPI_CHAR, target_cells[1], target_cells[1], grid_comm, &bot_recv);
    }
    
    // Left
    if (0 <= target_cells[2]) {
        get_cols(grid, out_buf_left, 0);
        MPI_Isend(out_buf_left, EDGE_LEN * loc_h, MPI_CHAR, target_cells[2], my_id, grid_comm, &left_send);
        MPI_Irecv(in_buf_left, EDGE_LEN * loc_h, MPI_CHAR, target_cells[2], target_cells[2], grid_comm, &left_recv);
    }
    
    // Right
    if (0 <= target_cells[3]) {
        get_cols(grid, out_buf_left, loc_w-EDGE_LEN);
        MPI_Isend(out_buf_left, EDGE_LEN * loc_h, MPI_CHAR, target_cells[3], my_id, grid_comm, &right_send);
        MPI_Irecv(in_buf_left, EDGE_LEN * loc_h, MPI_CHAR, target_cells[3], target_cells[3], grid_comm, &right_recv);
    }
    

    // Calculate all except borders
    #pragma omp parallel
    {
        #pragma omp for collapse(2) schedule(guided)
        for(int y = 1; y < loc_h-1; y++)
            for(int x = 1; x < loc_w-1; x++)
                next_grid[y][x] = next_state(grid,y,x);
    }

    // Wait for send to top
    if (0 <= target_cells[0]) {
        MPI_Wait(&top_send, MPI_STATUS_IGNORE);
    }
    // Wait for recv from bot and send to bot
    if (0 <= target_cells[1]) {
        MPI_Wait(&bot_recv, MPI_STATUS_IGNORE);
        set_rows(grid, in_buf_bot, loc_h);
        MPI_Wait(&bot_send, MPI_STATUS_IGNORE);
    }
    // Wait for recv from top
    if (0 <= target_cells[0]) {
        MPI_Wait(&top_recv, MPI_STATUS_IGNORE);
        set_rows(grid, in_buf_top, -EDGE_LEN);
    }

    // Wait for send to left
    if (0 <= target_cells[2]) {
        MPI_Wait(&left_send, MPI_STATUS_IGNORE);
    }
    // Wait for recv from right and send to right
    if (0 <= target_cells[3]) {
        MPI_Wait(&right_recv, MPI_STATUS_IGNORE);
        set_cols(grid, in_buf_right, loc_h);
        MPI_Wait(&right_send, MPI_STATUS_IGNORE);
    }// Wait for recv from left
    if (0 <= target_cells[2]) {
        MPI_Wait(&left_recv, MPI_STATUS_IGNORE);
        set_cols(grid, in_buf_left, -EDGE_LEN);
    }
    
    // Calculate all borders and update grid
    #pragma omp parallel
    {
        #pragma omp for collapse(2) schedule(guided)
        for(int x = -EDGE_LEN+1; x < loc_w+EDGE_LEN; x++) {
            for(int count = 0; count < EDGE_LEN; count++) {
                next_grid[0-count][x] = next_state(grid,0-count,x);
                next_grid[loc_h-1+count][x] = next_state(grid,loc_h-1+count,x);
            }
        }
        #pragma omp for collapse(2) schedule(guided)
        for(int y = 0; y < loc_h; y++) {
            for(int count = 0; count < EDGE_LEN; count++) {
                next_grid[y][0-count] = next_state(grid,y,0-count);
                next_grid[y][loc_h-1+count] = next_state(grid,y,loc_h-1+count);
            }
        }
        #pragma omp for collapse(2) schedule(guided)
        for(int y = 0; y < loc_h; y++)
            for(int x = 0; x < loc_w; x++)
                grid[y][x] = next_grid[y][x];
    }

    MPI_Barrier(grid_comm);
    free(out_buf_top);
    free(out_buf_bot);
    free(out_buf_left);
    free(out_buf_right);
    free(in_buf_top);
    free(in_buf_bot);
    free(in_buf_left);
    free(in_buf_right);
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
    int wrap_around[2] = {0,0};
    int dims;
    MPI_Comm grid_comm;

    int my_grid_id;
    int my_grid_coords[2];
    int my_neigbours[4];

    // ================== MPI ==================
    // Configure MPI parallelization based on input arguments
    MPI_Init(&argc, &argv);    
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
    // Set grid dimensions
    dims = set_dims(num_procs, dim_sizes);

    // Create communication for Cartesian topology
    MPI_Cart_create(MPI_COMM_WORLD, dims, dim_sizes, wrap_around, 0, &grid_comm);

    //Get coordinates
    MPI_Comm_rank(grid_comm, &my_grid_id);
    MPI_Cart_get(grid_comm, dims, dim_sizes, wrap_around, my_grid_coords);

    //Get neighbouring processes
    MPI_Cart_shift(grid_comm, 0, 1, &my_neigbours[2], &my_neigbours[3]);
    MPI_Cart_shift(grid_comm, 1, 1, &my_neigbours[0], &my_neigbours[1]);

    // Each process generates its own grid
    int h_size = (int)round((double) H / (double) dim_sizes[0]);
    int w_size = (int)round((double) W / (double) dim_sizes[1]);
    if (my_id == num_procs-1) {
        h_size = H - h_size * (dim_sizes[0]-1); //last process gets less
        w_size = W - w_size * (dim_sizes[1]-1);
    }

    if (my_grid_id == 0) 
        printf("Res:%dx%d | grid: %dx%d | edge: %d | iter: %d\n", W, H, dim_sizes[0], dim_sizes[1], EDGE_LEN, ITER);
    

    // Initialize local grids
    char** grid = alloc_grid(h_size, w_size, EDGE_LEN, dims == 2 ? EDGE_LEN : 0);
    char** next_grid = alloc_grid(h_size, w_size, EDGE_LEN, dims == 2 ? EDGE_LEN : 0);
    
    // Setup local sizes for process
    loc_h = h_size;
    loc_w = w_size;    
    proc_seed = my_id;
    printf("id:%d [%d,%d] | neighb: [%d,%d,%d,%d] | sizes: [%d,%d]\n", my_grid_id, my_grid_coords[0], my_grid_coords[1], my_neigbours[0], my_neigbours[1], my_neigbours[2], my_neigbours[3], loc_w, loc_h);

    init(grid, 0.7);
    if (my_id == 0) spark(grid, 10);

    double t_start = omp_get_wtime();
    
    // Simulation loop
    for (t = 0; t < ITER; t++) {
        // printf("id=%d\n", my_grid_id);
        // print_grid(grid, NULL);
        // printf("\n");
        if (t % EDGE_LEN == 0) {
            // Perform iteration with send
            if (dims == 2) 
                update_send_blocks(grid, next_grid, my_neigbours, grid_comm, my_grid_id);
            else 
                update_send(grid, next_grid, my_neigbours, grid_comm, my_grid_id);
        } else {
            // Perform iteration without send
            update_nosend(grid, next_grid, t % EDGE_LEN, dims == 2);
        }
    }
    double t_end = omp_get_wtime();

    if (my_grid_id == 0)
        printf("Elapsed time = %f\n", t_end-t_start);

    MPI_Finalize();
}


// Compile with: mpicc -lm -fopenmp parallel.c -o out
// Run with: sbatch --reservation=fri ./job.sh