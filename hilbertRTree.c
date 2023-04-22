#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M 4
#define m 2
#define MAX_CHILDREN M //M = 4; MAXIMUM NUMBER OF CHILDREN
#define MAX_MBR 2*M-1

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

typedef struct Node* NODE;
struct Node {
    int is_leaf;          // 1 if the node is a leaf, 0 otherwise
    int count;            // number of rectangles in the node
    Rectangle mbr[MAX_MBR];    // minimum bounding rectangles of the rectangles in the node
    float lhv;            // Hilbert value of the rectangles in the node
    struct Node *parent;  // pointer to the parent node
    struct Node *children[MAX_CHILDREN]; // pointers to the child nodes
};
