#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>


#define ITER 100//number of iterations
#define P_ST 0.01 //prob. of fire starting in a cell at init
#define P_BD 0.6 //prob. of burning cell burning down
#define P_H 0.58 //constant spreading prob.
#define IMAGE_PATH "./temp/" //image sequence is stored here

int H = 2000;       // height of grid
int W = 2000;       // width of grid
char symbols[5] = {'E',' ','.','&','*'};
int colors[5][3] = {{0,0,0},
                {80, 80, 80},
                {10, 40, 15},
                {255,215,80},
                {30, 30, 30 }};

/*  
    Cell states:
    0 - invalid state
    1 - non-flammable cell
    2 - flammable but intact
    3 - currently burning
    4 - had burned down
*/

// ================ UTILITY ================
// Get random number
double rnd(){
    return (double) rand() / (double) RAND_MAX;
}


// ================ GRID MANIPULATION ================ 
// Print the grid in current state
void print_grid(char** grid){
    for(int y = 0; y < H; y++){
        for(int x = 0; x < W; x++){
            printf("%c", symbols[grid[y][x]]);
        }
        printf("\n");
    }         
}

// Returns grid that can be accessed like [i][j], and has padded edges allocated
char** alloc_grid(int h, int w){
    char* alloc = (char *) malloc(h * w * sizeof(char));
    char** grid = (char **) malloc(h * sizeof(char *));
    for(int i = 0; i < h; i++)
            grid[i] = &alloc[i*w];
    return grid;
}



// ================ CELL SIMULATION ================
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
        x0 = rand() % W;
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
                if (y+i < 0 || y+i >= H || x+j < 0 || x+j >= W)
                    continue;
                if (grid[y+i][x+j] == 3 && rnd() < p_spread)
                    return 3;
            }

    }
    return grid[y][x];
}



// ================ UPDATE METHODS ================
// Update grid
void update(char** grid, char** next_grid){ 
    for(int y = 0; y < H; y++)
        for(int x = 0; x < W; x++)
            next_grid[y][x] = next_state(grid, y,x);

    for(int y = 0; y < H; y++)
        for(int x = 0; x < W; x++)
            grid[y][x] = next_grid[y][x];
}

int main(int argc, char* argv[]){

    if (argc <= 2) {
        printf("Argument order: H W EDGE_LEN BLOCKS\n");
        exit(1);
    }

    H = atoi(argv[1]);
    W = atoi(argv[2]);

    // Check if all arguments are valid
    if (H < 32 || 8000 < H) {
        printf("H:%d too small or too big", H);
        exit(1);
    }
    if (W < 32 || 8000 < W) {
        printf("W:%d too small or too big\n", W);
        exit(1);
    }
    

    // Initialize local grids
    char** grid = alloc_grid(H, W);
    char** next_grid = alloc_grid(H, W);
    
    init(grid, 0.7);
    spark(grid, 5);

    clock_t t_start = clock();
    // Simulation loop
    for (int t = 0; t < ITER; t++) {
        update(grid, next_grid);
    }
    clock_t t_end = clock();
    printf("%f\n", (float)(t_end - t_start) / CLOCKS_PER_SEC);
}
