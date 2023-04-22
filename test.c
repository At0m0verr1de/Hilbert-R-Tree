#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

// EVERY NODE HAS BETWEEN m and M entries and children unless it is root
#define M 4
#define m 2
#define MAX_CHILDREN M //M = 4; MAXIMUM NUMBER OF CHILDREN
#define MIN_CHILDREN m
#define MAX_POINTS 21
int CURRENT_ID = 1;
int num_results = 0;

typedef struct Point Point;
struct Point{
    double x,y;
};
typedef struct Rectangle Rectangle;
struct Rectangle {
    Point bottom_left;
    Point top_right; // coordinates of the rectangle
    float h;              // Hilbert value of the rectangle center
};
typedef struct LeafEntry LeafEntry;
struct LeafEntry {
    Rectangle mbr;  // Minimum bounding rectangle
    int obj_id;  // Pointer to the object description record
};

struct LeafNode {
    int num_entries;
    struct LeafEntry entries[MAX_CHILDREN];
};

typedef struct NonLeafEntry NonLeafEntry;
struct NonLeafEntry {
    Rectangle mbr;  // Minimum bounding rectangle
    struct Node* child_ptr;  // Pointer to the child node
    int largest_hilbert_value;  // Largest hilbert value among the data records enclosed by the MBR
};

struct NonLeafNode {
    int num_entries;
    struct NonLeafEntry entries[MAX_CHILDREN];
    
};
typedef struct Node* NODE;
struct Node {
    int is_leaf;
    NODE parent_ptr;
    union {
        struct NonLeafNode non_leaf_node;  // Non-leaf node
        struct LeafNode leaf_node;  // Leaf node
    } u;
};

typedef struct HilbertRTree hilbertRTree;
struct HilbertRTree {
    int height;  // Height of the tree
    struct Node* root;  // Root node of the tree
};
void rotate2D(int *x, int *y, int rx, int ry, int order) {
    if (!ry) {
        if (rx) {
            *x = order-1 - *x;
            *y = order-1 - *y;
        }
        int t = *x;
        *x = *y;
        *y = t;
    }
}

int hilbert_2d(int x, int y, int order) {
    int hilbert_value = 0;
    int mask = (1 << order) - 1;

    for (int s = order - 1; s >= 0; s--) {
        int rx = (x >> s) & 1;
        int ry = (y >> s) & 1;
        hilbert_value += mask * mask * ((3 * rx) ^ ry);
        rotate2D(&x, &y, rx, ry, 1 << s);
        mask >>= 1;
    }

    return hilbert_value;
}


int main() {
   int x = 43;
int y = 70;
int order = 2;

int hilbert_value = hilbert_2d(x, y, order);
printf("The Hilbert value of (%d, %d) with order %d is %d\n", x, y, order, hilbert_value);

}