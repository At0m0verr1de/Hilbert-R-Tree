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
    int x,y;
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
    NODE parent_ptr;
};

struct NonLeafEntry {
    Rectangle mbr;  // Minimum bounding rectangle
    struct Node* child_ptr;  // Pointer to the child node
    int largest_hilbert_value;  // Largest hilbert value among the data records enclosed by the MBR
};

struct NonLeafNode {
    int num_entries;
    struct NonLeafEntry entries[MAX_CHILDREN];
    NODE parent_ptr;
};
typedef struct Node* NODE;
struct Node {
    int is_leaf;
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

//IS_LEAF = 1 IF LEAF NODE. M = MAXIMUM NUMBER OF CHILDREN THAT NODE CAN HAVE

// Function to calculate the minimum bounding rectangle that contains two given rectangles
Rectangle calculateMBR(Rectangle r1, Rectangle r2){
    Point bottom_leftNew;
    Point top_rightNew;
    if(r1.bottom_left.x <= r2.bottom_left.x){
        bottom_leftNew.x = r1.bottom_left.x;
    }
    else{
        bottom_leftNew.x = r2.bottom_left.x;
    }
    if(r1.bottom_left.y <= r2.bottom_left.y){
        bottom_leftNew.y = r1.bottom_left.y;
    }
    else{
        bottom_leftNew.y = r2.bottom_left.y;
    }

    if(r1.top_right.x <= r2.top_right.x){
        top_rightNew.x = r2.top_right.x;
    }
    else{
        top_rightNew.x = r1.top_right.x;
    }
    if(r1.top_right.y <= r2.top_right.y){
        top_rightNew.y = r2.top_right.y;
    }
    else{
        top_rightNew.y = r1.top_right.y;
    }

    Rectangle new_rect;
    new_rect.bottom_left = bottom_leftNew;
    new_rect.top_right = top_rightNew;
    return new_rect;
}


// PRINT THE MBR - TOP RIGHT POINT AND BOTTOM LEFT POINT
void printMBR(Rectangle rect){
    printf("MBR = %d %d %d %d\n", rect.bottom_left.x, rect.bottom_left.y, rect.top_right.x, rect.top_right.y);
    return;
}

// Compute the Hilbert value for a given rectangle defined by its bottom-left and top-right corners
unsigned long long HilbertValue(Point bottom_left, Point top_right) {
    // Determine the length of the longest edge of the rectangle
    double max_side = fmax(top_right.x - bottom_left.x, top_right.y - bottom_left.y);
    
    // Compute the number of levels in the Hilbert curve
    int num_levels = ceil(log2(max_side));
    
    // Compute the number of cells in the Hilbert curve
    unsigned long long num_cells = pow(2, 2*num_levels);
    
    // Compute the Hilbert value for the given rectangle
    unsigned long long hilbert_value = 0;
    for (int i = num_levels - 1; i >= 0; i--) {
        unsigned long long mask = 1 << i;
        double x = (bottom_left.x & mask) >> i;
        double y = (bottom_left.y & mask) >> i;
        hilbert_value += mask * pow((3 * x), y);
        
        if (y == 0) {
            if (x == 1) {
                bottom_left.x = max_side - bottom_left.x;
                top_right.x = max_side - top_right.x;
            }
            double tmp = bottom_left.x;
            bottom_left.x = bottom_left.y;
            bottom_left.y = tmp;
            tmp = top_right.x;
            top_right.x = top_right.y;
            top_right.y = tmp;
        }
    }
    
    return hilbert_value;
}


NODE ChooseLeaf(NODE n, Rectangle r, int h){
    /* RETURNS THE LEAF NODE IN WHICH TO PLACE A NEW RECTANGE*/
    
    /* IF N IS A LEAF, RETURN N*/
    if(n->is_leaf == 1){
        return n;
    }

    /* IF N IS A NON-LEAF NODE, CHOOSE ENTRY (R, PTR, LHV) WITH MINIMUM LHV
    GREATER THAN H*/
    float min_LHV= INFINITY;
    NODE next_node = NULL;
    for(int i = 0; i<n->u.non_leaf_node.num_entries; i++){
            if(n->u.non_leaf_node.entries[i].largest_hilbert_value > h && n->u.non_leaf_node.entries[i].largest_hilbert_value < min_LHV){
                min_LHV = n->u.non_leaf_node.entries[i].largest_hilbert_value;
                next_node = n->u.non_leaf_node.entries[i].child_ptr;
            }
        
    }
    // IF ALL CHILDREN HAVE LHV LESS THAN H
    if(next_node == NULL){
        //CHOOSE THE CHILD NODE WITH LARGEST LHV
        for(int i = 0; i<n->u.non_leaf_node.num_entries; i++){
            
                if(n->u.non_leaf_node.entries[i].largest_hilbert_value > min_LHV){
                    min_LHV = n->u.non_leaf_node.entries[i].largest_hilbert_value;
                    next_node = n->u.non_leaf_node.entries[i].child_ptr;
                }
            
        }
    }

    /* DESCEND UNTIL A LEAF NODE IS REACHED*/
    return ChooseLeaf(next_node, r, h);
}

void Insert(NODE root, Rectangle rectangle){
    NODE leafNode = ChooseLeaf(root, rectangle, rectangle.h);
    NODE newLeafNode = NULL;
    if(leafNode->u.leaf_node.num_entries < M){
        //LEAF NODE HAS SPACE
        int i = 0;
        while(i < leafNode->u.leaf_node.num_entries && leafNode->u.leaf_node.entries[i].mbr.h < rectangle.h){
            i++;
        }

        for(int j = leafNode->u.leaf_node.num_entries; j>i; j--){
            leafNode->u.leaf_node.entries[j] = leafNode->u.leaf_node.entries[j-1];
        }
        LeafEntry entry;
        entry.mbr = rectangle; 
        entry.obj_id = ++CURRENT_ID;
        leafNode->u.leaf_node.entries[i] = entry;
        leafNode->u.leaf_node.num_entries++;
    }
    else{
        //LEAF NODE IS FULL
        newLeafNode = HandleOverFlow(leafNode, rectangle);
        //RETURNS THE NEW LEAF IF SPLIT WAS INEVITABLE
    }

    //Propogate changes upward

    //FORM A SET S CONTAINING L: COOPERATING SIBLINGS AND NEW LEAF (IF ANY)
     NODE* S = (NODE*)malloc(sizeof(NODE) * MAX_POINTS);
    for (int i = 0; i < MAX_POINTS; i++) {
        S[i] = NULL;
    }
    S[0] = leafNode;
    if (newLeafNode) {
        S[1] = newLeafNode;
    }
    int numSiblings = 1;
    NODE parentNode = leafNode->u.non_leaf_node.parent_ptr;
    //GO FROM CURRENT NODE TO ROOT FINDING SIBLINGS
    //WITH LESS THAN MAXIMUM POINTERS
    while (parentNode) {
            int index = -1;
            for (int i = 0; i < parentNode->u.non_leaf_node.num_entries; i++) {
                if (parentNode->u.non_leaf_node.entries[i].child_ptr == S[0]) {
                    index = i;
                    break;
                }
            }
            if (index > 0 && parentNode->u.non_leaf_node.entries[index - 1].child_ptr->u.leaf_node.num_entries < M) {
                S[++numSiblings] = parentNode->u.non_leaf_node.entries[index - 1].child_ptr;
            } else if (index == 0 && parentNode->u.non_leaf_node.entries[index + 1].child_ptr->u.leaf_node.num_entries < M) {
                S[++numSiblings] = parentNode->u.non_leaf_node.entries[index + 1].child_ptr;
            }
            S[0] = parentNode;
            parentNode = parentNode->u.non_leaf_node.parent_ptr;
    }
    //HOW TO ADD COOPERATING SIBLINGS?
    AdjustTree(S);

    //IF NODE SPLIT CAUSED ROOT TO SPLIT, CREATEA NEWROOT WITH CHILDREN
    //AS RESULTING NODES
    // Check if the root split
if (S[0]->parent_ptr == NULL) {
    NODE newRoot = (NODE) malloc(sizeof(struct Node));
    newRoot->is_leaf = 0;
    newRoot->u.non_leaf_node.num_entries = 1;
    newRoot->u.non_leaf_node.parent_ptr = NULL;
    
    newRoot->u.non_leaf_node.entries[0].child_ptr = S[0];
    newRoot->u.non_leaf_node.entries[0].mbr = (S[0]->u.);
    
    newRoot->u.non_leaf_node.entries[1].child_ptr = S[1];
    newRoot->u.non_leaf_node.entries[1].mbr = CalculateMBR(S[1]);

    S[0]->parent_ptr = newRoot;
    S[1]->parent_ptr = newRoot;

    root = newRoot;
}

}


int main() {
    FILE *fp;
    Point points[MAX_POINTS];
    int num_points = 0;

    // Open the file containing the data points
    fp = fopen("data.txt", "r");
    if (fp == NULL) {
        printf("Error opening file.\n");
        return 1;
    }

    // Read the data points from the file
    while (fscanf(fp, "%d %d\n", &points[num_points].x, &points[num_points].y) == 2) {
        num_points++;
        printf("POINT  = %d %d\n", points[num_points-1].x, points[num_points-1].y);
    }

    // Close the file
    fclose(fp);

    //LETS ASSUME THE RECTANGLES INITIALLY ARE THE DATA POINTS
    //AND THEIR CENTRES ARE THE POINTS TOO

    // CREATE RECTANGLES
    Rectangle rectangles[MAX_POINTS];
    for(int i = 0; i<num_points; i++){
        rectangles[i].bottom_left = points[i];
        rectangles[i].top_right = points[i];
        rectangles[i].h = hilbert_value(points[i]);        
    }


    
    NODE root = (NODE)malloc(sizeof(struct Node));
    root->is_leaf = 1; //INITIALLY ROOT IS A LEAF NODE

    //CREATE LEAF ENTRIES AND INSERT TO ROOT NODE
    for(int i = 0; i<MAX_POINTS; i++){
        Insert(root, rectangles[i]); //INSERT(NODE ROOT, RECTANGLE R)
    }
    

    return 0;
}


// SEARCH ALGORITHM
//-> NONLEAF - THOSE WITH MBR INTERSECTING THE QUERY WINDOW W
//-> LEAF - THOSE WITH MBR INTERSECTING THE QUERY WINDOW W
int num_results = 0;void search(NODE root, Rectangle rectangle, LeafEntry* results){
    if(root->is_leaf == 1){
        for(int i = 0 ; i<root->u.leaf_node.num_entries; i++){
            if(intersects(root->u.leaf_node.entries[i].mbr, rectangle)){
                results[num_results++] = root->u.leaf_node.entries[i];
            }
        }
    }
    else{
        for(int i = 0; i<root->u.non_leaf_node.num_entries; i++){
            if(intersects(root->u.non_leaf_node.entries[i].mbr, rectangle)){
                search(root->u.non_leaf_node.entries[i].child_ptr, rectangle, results);
            }
        }
    }
}

bool intersects(Rectangle r1, Rectangle r2){
    return !(
        r1.top_right.x < r2.bottom_left.x || r2.top_right.x < r1.bottom_left.x || 
        r1.top_right.y < r2.bottom_left.y || r2.top_right.y < r1.bottom_left.y
    );
}
