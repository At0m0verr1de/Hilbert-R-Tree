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


// Helper function to calculate the Hilbert value of a point
uint32_t interleave(uint32_t x)
{
    x = (x | (x << 8)) & 0x00FF00FF;
    x = (x | (x << 4)) & 0x0F0F0F0F;
    x = (x | (x << 2)) & 0x33333333;
    x = (x | (x << 1)) & 0x55555555;
    return x;
}

// Calculate the Hilbert value of a point
uint32_t hilbertXYToIndex(uint32_t n, uint32_t x, uint32_t y)
{
    x = x << (16 - n);
    y = y << (16 - n);

    uint32_t A, B, C, D;

    // Initial prefix scan round, prime with x and y
    {
        uint32_t a = x ^ y;
        uint32_t b = 0xFFFF ^ a;
        uint32_t c = 0xFFFF ^ (x | y);
        uint32_t d = x & (y ^ 0xFFFF);

        A = a | (b >> 1);
        B = (a >> 1) ^ a;

        C = ((c >> 1) ^ (b & (d >> 1))) ^ c;
        D = ((a & (c >> 1)) ^ (d >> 1)) ^ d;
    }

    {
        uint32_t a = A;
        uint32_t b = B;
        uint32_t c = C;
        uint32_t d = D;

        A = ((a & (a >> 2)) ^ (b & (b >> 2)));
        B = ((a & (b >> 2)) ^ (b & ((a ^ b) >> 2)));

        C ^= ((a & (c >> 2)) ^ (b & (d >> 2)));
        D ^= ((b & (c >> 2)) ^ ((a ^ b) & (d >> 2)));
    }

    {
        uint32_t a = A;
        uint32_t b = B;
        uint32_t c = C;
        uint32_t d = D;

        A = ((a & (a >> 4)) ^ (b & (b >> 4)));
        B = ((a & (b >> 4)) ^ (b & ((a ^ b) >> 4)));

        C ^= ((a & (c >> 4)) ^ (b & (d >> 4)));
        D ^= ((b & (c >> 4)) ^ ((a ^ b) & (d >> 4)));
    }

    // Final round and projection
    {
        uint32_t a = A;
        uint32_t b = B;
        uint32_t c = C;
        uint32_t d = D;

        C ^= ((a & (c >> 8)) ^ (b & (d >> 8)));
        D ^= ((b & (c >> 8)) ^ ((a ^ b) & (d >> 8)));
    }

    // Undo transformation prefix scan
    uint32_t a = C ^ (C >> 1);
    uint32_t b = D ^ (D >> 1);

    // Recover index bits
    uint32_t i0 = x ^ y;
    uint32_t i1 = b | (0xFFFF ^ (i0 | a));

    return ((interleave(i1) << 1) | interleave(i0)) >> (32 - 2 * n);
}

// Calculate the Hilbert value of a rectangle
uint32_t hilbertValue(Rectangle rect){
    uint32_t x = (rect.bottom_left.x + rect.top_right.x)/2;
    uint32_t y = (rect.bottom_left.y + rect.top_right.y)/2;
    uint32_t hilbert_value = hilbertXYToIndex(16, x, y);
    return hilbert_value;
}

