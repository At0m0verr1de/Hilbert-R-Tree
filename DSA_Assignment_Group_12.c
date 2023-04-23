#include <stdint.h>
#include <stdio.h>

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

typedef struct HilbertRTree HilbertRTree;
struct HilbertRTree {
    int height;  // Height of the tree
    struct Node* root;  // Root node of the tree
};

//FUNCTIONS TO CREATE NODE AND TREE
HilbertRTree* new_hilbertRTree(){
    HilbertRTree* tree = (HilbertRTree*)malloc(sizeof(HilbertRTree));
    tree->root = NULL;
    tree->height = 0;
    return tree;
}

NODE new_node(int is_leaf){
    NODE node = (NODE)malloc(sizeof(struct Node));
    node->is_leaf = is_leaf;
    node->parent_ptr = NULL;
    if(is_leaf){
        node->u.leaf_node.num_entries = 0;
    }
    else{
        node->u.non_leaf_node.num_entries = 0;
    }
    return node;
}


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


// Print the MBR - Top Right and Bottom Left coordinates
void printMBR(Rectangle rect){
    printf("MBR = %d %d %d %d\n", rect.bottom_left.x, rect.bottom_left.y, rect.top_right.x, rect.top_right.y);
    return;
}


