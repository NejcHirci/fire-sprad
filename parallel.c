#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "FreeImage.h"
#include "omp.h"
#include "mpi.h"

#define H 64
#define W 64
#define ITER 64 //number of iterations
#define P_ST 0.01 //prob. of fire starting in a cell at init
#define P_BD 0.6 //prob. of burning cell burning down
#define P_H 0.58 //constant spreading prob.
#define IMAGE_PATH "./temp/" //image sequence is stored here

int grid[H][W];
int next_grid[H][W];
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
void init(char* grid){
    for(int y = 0; y < H; y++)
        for(int x = 0; x < W; x++){
            grid[y*H+x] = rand() % 2;
        }
}

// Start the fire in some cells: n - number of sparks
void spark(char* grid, int n){
    int numSparks = 0;
    int y0, x0;
    while (numSparks < n){
        y0 = rand() % H;
        x0 = rand() % W;
        grid[y0*H+W] = 3;
    }
}

// GET STATE OF A CELL AT t+1
int next_state(int y, int x, int t){

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
void update(){
    for(int y = 0; y < H; y++)
        for(int x = 0; x < W; x++)
            next_grid[y][x] = next_state(x,y,t);

    for(int y = 0; y < H; y++)
        for(int x = 0; x < W; x++)
            grid[y][x] = next_grid[y][x];
}

// Parallel state update for whole grid
void updateParallel(){
    #pragma omp parallel for schedule(guided) private(sum, sum2) num_threads(NTHREADS)
    for(int y = 0; y < H; y++)
        for(int x = 0; x < W; x++)
            next_grid[y][x] = next_state(x,y,t);
    
}

// Print the grid in current state
void print_grid(){
    for(int y = 0; y < H; y++){
        for(int x = 0; x < W; x++){
            printf("%c", symbols[grid[y][x]]);
        }
        printf("\n");
    }
}

// Render grid as png image
void render_grid(){
    int pitch = ((32 * W + 31) / 32) * 4;

    for(int y = 0; y < H; y++)
        for(int x = 0; x < W; x++){
            int state = grid[y][x];
            image[4*y*W + 4*x + 0] = colors[state][2];  // Blue
            image[4*y*W + 4*x + 1] = colors[state][1];  // Green
			image[4*y*W + 4*x + 2] = colors[state][0];  // Red
			image[4*y*W + 4*x + 3] = 255;               // Alpha
        }

    char image_name[16];
    sprintf(image_name, "%s%d.png", IMAGE_PATH, t);

    FIBITMAP *dst = FreeImage_ConvertFromRawBits(image, W, H, pitch,
		32, FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK, TRUE);
	FreeImage_Save(FIF_PNG, dst, image_name, 0);
}

int main(int argc, char* argv[]){

    int myid, numOfProc
    MPI_Status status;

    char* globalGrid = (char *)malloc(H * W * sizeof(unsigned char));

    // Initialize grid and spark fire in some cells
    init(globalGrid);
    spark(globalGrid, 20);

    // ================== MPI ==================

    // Configure MPI paralellization based on input arguments
    MPI_Init(&argc, &argv);    
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);

}