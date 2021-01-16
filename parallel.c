#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "FreeImage.h"
#include "omp.h"
//#include "mpi.h"

#define H 32
#define W 32
#define ITER 20//number of iterations
#define P_ST 0.01 //prob. of fire starting in a cell at init
#define P_BD 0.6 //prob. of burning cell burning down
#define P_H 0.58 //constant spreading prob.
#define IMAGE_PATH "./temp/" //image sequence is stored here

#define X_EDGE 1 //width of the horizontal "edge" of neighbouring cells to exhange
#define Y_EDGE 0 //width of the vertical edge (has to be 0 or same as X_edge)

//int grid[H][W];
//int next_grid[H][W];
int t;

char symbols[5] = {' ',' ','.','&','*'};
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
    return (double) rand() / (double) RAND_MAX;
}

// Initialize cells (to 1 or 2)
void init(char** grid, double density){
    for(int y = 0; y < H; y++)
        for(int x = 0; x < W; x++){
            grid[y][x] = rnd() < density ? 2 : 1;
        }
}

// Start the fire in some cells: n - number of sparks
void spark(char** grid, int n){
    int y0, x0;
    for(int i = 0; i < n; i++){
        y0 = rand() % H;
        x0 = rand() % W ;
        grid[y0][x0] = 3;
    }
}

// GET STATE OF A CELL AT t+1
int next_state(char** grid, int y, int x, int t){

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
                if (y+i < 0 || y+i >= H || x+j < 0 || x+j >= W)
                    continue;
                if (grid[y+i][x+j] == 3 && rnd() < p_spread)
                    return 3;
            }

    }

    return grid[x][y];
}

// Update states for whole grid (one iteration)
// void update(){
//     for(int y = 0; y < H; y++)
//         for(int x = 0; x < W; x++)
//             next_grid[y][x] = next_state(x,y,t);

//     for(int y = 0; y < H; y++)
//         for(int x = 0; x < W; x++)
//             grid[y][x] = next_grid[y][x];
// }

// Parallel state update for whole grid
void updateParallel(char** grid, char** next_grid, int li){
    int e = X_EDGE - 1 - li;

    // najhitreje je z -c16
    #pragma omp parallel
    {

        #pragma omp for collapse(2) schedule(static)
        for(int y = -e; y < H+e; y++)
            for(int x = -e; x < W+e; x++)
                next_grid[y][x] = next_state(grid,x,y,t);

        #pragma omp for collapse(2) schedule(static)
        for(int y = -e; y < H+e; y++)
            for(int x = -e; x < W+e; x++)
                grid[y][x] = next_grid[y][x];
    }
}


// Print the grid in current state
void print_grid(char** grid){
    for(int y = 0; y < H; y++){
        for(int x = 0; x < W; x++){
            printf("%c", symbols[grid[y][x]]);
        }
        printf("\n");
    }
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

// Exhange edge data with parallel processes
void exchange_data(char** grid){

}

void send_xedge(){
    //grid edge to buffer and send
}

void send_yedge(){
    //grid edge to buffer and send
}

void recv_xedge(){
    //receive buffer, buffer to grid edge
}

void recv_yedge(){
    //receive buffer, buffer to grid edge
}

void send_corner(){

}

void recv_corner(){

}


int main(int argc, char* argv[]){

    int P, myid;
    //MPI_Status status;

    // ! has to be allocated in process #1 and saved to a shared file or smh
    //char* globalGrid = (char *)malloc(H * W * sizeof(unsigned char));

    // Allocate the "band" for this process - sneaky version - upam a OMP parallel dela s tem
    // char* local_grid = (char *) malloc(H * W * sizeof(char))


    char** grid = alloc_grid(H, W, X_EDGE, Y_EDGE);
    char** next_grid = alloc_grid(H, W, X_EDGE, Y_EDGE);

    // Initialize grid and spark fire in some cells
    init(grid, 0.7);
    spark(grid, 3);


    // ================== MPI ==================

    // Configure MPI parallelization based on input arguments
    // MPI_Init(&argc, &argv);    
    // MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    // MPI_Comm_size(MPI_COMM_WORLD, &P);

    
    // 1. Exchange edge data
    // 2. Run local iterations
    // 3. Repeat

    double t_start = omp_get_wtime();

    for(t = 0; t < ITER; t++){
        if (t % X_EDGE == 0){

            // EXHANGE DATA WITH MPI
        }

        print_grid(grid); 
        printf("\n");
        //render_grid(); 

        int li = t % X_EDGE; //local iterations since the last edge exchange
        updateParallel(grid, next_grid, li);
    }

    double t_end = omp_get_wtime();
    printf("Elapsed time = %f\n", t_end-t_start);
}