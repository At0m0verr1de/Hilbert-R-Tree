#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
// #include <hilbertcurve/hilbertcurve.h>
// EVERY NODE HAS BETWEEN m and M entries and children unless it is root
#define M 4
#define m 2
#define MAX_CHILDREN M // M = 4; MAXIMUM NUMBER OF CHILDREN
#define MIN_CHILDREN m
#define MAX_POINTS 21

typedef unsigned long long ull;
int CURRENT_ID = 1;
int num_results = 0;

typedef struct Point Point;
struct Point
{
    int x, y;
};

typedef struct Rectangle Rectangle;
struct Rectangle
{
    Point bottom_left;
    Point top_right; // coordinates of the rectangle
    float h;         // Hilbert value of the rectangle center
};

typedef struct LeafEntry LeafEntry;
struct LeafEntry
{
    Rectangle mbr; // Minimum bounding rectangle
    int obj_id;    // Pointer to the object description record
};

struct LeafNode
{
    int num_entries;
    struct LeafEntry entries[MAX_CHILDREN];
};

typedef struct NonLeafEntry NonLeafEntry;
struct NonLeafEntry
{
    Rectangle mbr;          // Minimum bounding rectangle
    struct Node *child_ptr; // Pointer to the child node
    int lhv;                // Largest hilbert value among the data records enclosed by the MBR
};

struct NonLeafNode
{
    int num_entries;
    struct NonLeafEntry entries[MAX_CHILDREN];
};
typedef struct Node *NODE;
struct Node
{
    int is_leaf;
    NODE parent_ptr;
    union
    {
        struct NonLeafNode non_leaf_node; // Non-leaf node
        struct LeafNode leaf_node;        // Leaf node
    } u;
};

typedef struct HilbertRTree HilbertRTree;
struct HilbertRTree
{
    int height;        // Height of the tree
    struct Node *root; // Root node of the tree
};

// FUNCTIONS TO CREATE NODE AND TREE
HilbertRTree *new_hilbertRTree()
{
    HilbertRTree *tree = (HilbertRTree *)malloc(sizeof(HilbertRTree));
    tree->root = NULL;
    tree->height = 0;
    return tree;
}

NODE new_node(int is_leaf)
{
    NODE node = (NODE)malloc(sizeof(struct Node));
    node->is_leaf = is_leaf;
    node->parent_ptr = NULL;
    node->u.leaf_node.num_entries = 0;
    node->u.non_leaf_node.num_entries = 0;
    return node;
}

// IS_LEAF = 1 IF LEAF NODE. M = MAXIMUM NUMBER OF CHILDREN THAT NODE CAN HAVE

// Function to calculate the minimum bounding rectangle that contains two given rectangles
Rectangle calculateMBR(Rectangle r1, Rectangle r2)
{
    Point bottom_leftNew;
    Point top_rightNew;
    if (r1.bottom_left.x <= r2.bottom_left.x)
    {
        bottom_leftNew.x = r1.bottom_left.x;
    }
    else
    {
        bottom_leftNew.x = r2.bottom_left.x;
    }
    if (r1.bottom_left.y <= r2.bottom_left.y)
    {
        bottom_leftNew.y = r1.bottom_left.y;
    }
    else
    {
        bottom_leftNew.y = r2.bottom_left.y;
    }

    if (r1.top_right.x <= r2.top_right.x)
    {
        top_rightNew.x = r2.top_right.x;
    }
    else
    {
        top_rightNew.x = r1.top_right.x;
    }
    if (r1.top_right.y <= r2.top_right.y)
    {
        top_rightNew.y = r2.top_right.y;
    }
    else
    {
        top_rightNew.y = r1.top_right.y;
    }

    Rectangle new_rect;
    new_rect.bottom_left = bottom_leftNew;
    new_rect.top_right = top_rightNew;
    return new_rect;
}

// PRINT THE MBR - TOP RIGHT POINT AND BOTTOM LEFT POINT
void printMBR(Rectangle rect)
{
    printf("MBR = BLx:%d BLy:%d TRx:%d TRy:%d\n", rect.bottom_left.x, rect.bottom_left.y, rect.top_right.x, rect.top_right.y);
    return;
}

void InsertNode(NODE parent, NODE newNode)
{
    // newNode->is_leaf=0;
    // parent->is_leaf=0;
    newNode->parent_ptr = parent;                                                                                                                                       // SET PARENT POINTER OF NEW NODE TO PARENT
    parent->u.non_leaf_node.entries[parent->u.non_leaf_node.num_entries].child_ptr = newNode;                                                                           // SET CHILD POINTER OF PARENT TO NEW NODE
    parent->u.non_leaf_node.entries[parent->u.non_leaf_node.num_entries].mbr = calculateEntryMBR(parent->u.non_leaf_node.entries[parent->u.non_leaf_node.num_entries]); // SET MBR OF PARENT TO NEW NODE
    parent->u.non_leaf_node.entries[parent->u.non_leaf_node.num_entries].lhv = calculateLHV(parent->u.non_leaf_node.entries[parent->u.non_leaf_node.num_entries]);      // SET LHV OF PARENT TO NEW NODE
    parent->u.non_leaf_node.num_entries++;
}

// Compute the Hilbert value for a given rectangle defined by its bottom-left and top-right corners
ull HilbertValue(Point bottom_left, Point top_right)
{
    // Determine the length of the longest edge of the rectangle
    double max_side = fmax(top_right.x - bottom_left.x, top_right.y - bottom_left.y);

    // Compute the number of levels in the Hilbert curve
    int num_levels = ceil(log2(max_side));

    // Compute the number of cells in the Hilbert curve
    ull num_cells = pow(2, 2 * num_levels);

    // Compute the Hilbert value for the given rectangle
    ull hilbert_value = 0;
    for (int i = num_levels - 1; i >= 0; i--)
    {
        ull mask = 1 << i;
        double x = (bottom_left.x & mask) >> i;
        double y = (bottom_left.y & mask) >> i;
        hilbert_value += mask * pow((3 * x), y);

        if (y == 0)
        {
            if (x == 1)
            {
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
bool rectangles_equal(Rectangle *rect1, Rectangle *rect2)
{
    if (rect1->bottom_left.x != rect2->bottom_left.x ||
        rect1->bottom_left.y != rect2->bottom_left.y ||
        rect1->top_right.x != rect2->top_right.x ||
        rect1->top_right.y != rect2->top_right.y ||
        rect1->h != rect2->h)
    {
        return false;
    }
    return true;
}

bool nodes_equal(NODE node1, NODE node2)
{
    if (node1->is_leaf != node2->is_leaf)
    {
        return false;
    }
    if (node1->parent_ptr != node2->parent_ptr)
    {
        return false;
    }
    if (node1->is_leaf)
    {
        // Compare LeafNode fields
        struct LeafNode *leaf1 = &node1->u.leaf_node;
        struct LeafNode *leaf2 = &node2->u.leaf_node;
        if (leaf1->num_entries != leaf2->num_entries)
        {
            return false;
        }
        for (int i = 0; i < leaf1->num_entries; i++)
        {
            if (rectangles_equal(&leaf1->entries[i].mbr, &leaf2->entries[i].mbr) &&
                leaf1->entries[i].obj_id == leaf2->entries[i].obj_id)
            {
                continue;
            }
            else
            {
                return false;
            }
        }
    }
    else
    {
        // Compare NonLeafNode fields
        struct NonLeafNode *non_leaf1 = &node1->u.non_leaf_node;
        struct NonLeafNode *non_leaf2 = &node2->u.non_leaf_node;
        if (non_leaf1->num_entries != non_leaf2->num_entries)
        {
            return false;
        }
        for (int i = 0; i < non_leaf1->num_entries; i++)
        {
            if (rectangles_equal(&non_leaf1->entries[i].mbr, &non_leaf2->entries[i].mbr) &&
                non_leaf1->entries[i].child_ptr == non_leaf2->entries[i].child_ptr &&
                non_leaf1->entries[i].lhv == non_leaf2->entries[i].lhv)
            {
                continue;
            }
            else
            {
                return false;
            }
        }
    }
    return true;
}

// HELPER FUNCTION TO CHECK IN ARRAY
bool isInArray(NODE *arr, int size, NODE node)
{
    for (int i = 0; i < size; i++)
    {
        if (nodes_equal(arr[i], node))
        {
            return true;
        }
    }
    return false;
}
void rotate(uint32_t n, double *x, double *y, uint64_t rx, uint64_t ry)
{
    if (ry == 0)
    {
        if (rx == 1)
        {
            *x = n - 1 - *x;
            *y = n - 1 - *y;
        }
        double temp = *x;
        *x = *y;
        *y = temp;
    }
}

uint64_t xy2d(uint32_t n, double x, double y)
{
    uint64_t d = 0;
    for (int s = n / 2; s > 0; s /= 2)
    {
        uint64_t rx = (uint64_t)(x >= s);
        uint64_t ry = (uint64_t)(y >= s);
        d += s * s * ((3 * rx) ^ ry);
        rotate(s, &x, &y, rx, ry);
    }
    return d;
}

// Assumes the maximum x and y coordinates are 2^31-1
#define MAX_COORDINATE 2147483647
uint32_t get_hilbert_order(uint32_t num_dims, uint32_t bits_per_dim)
{
    uint32_t max_dim_index = (1 << bits_per_dim) - 1;
    uint32_t max_index = max_dim_index * num_dims;
    uint32_t order = 1;
    while (order < max_index)
    {
        order <<= 1;
    }
    return order;
}
uint64_t calculate_hilbert_value(double x, double y)
{
    // Scale the x and y coordinates to the range [0, 1]
    double scaled_x = x / MAX_COORDINATE;
    double scaled_y = y / MAX_COORDINATE;

    // Get the Hilbert curve order that fits the coordinates
    uint32_t order = get_hilbert_order(MAX_COORDINATE, 2);

    // Calculate the Hilbert value of the scaled point
    uint64_t hilbert_value = xy2d(order, scaled_x, scaled_y);

    return hilbert_value;
}

/*ADJUST TREE ASCEND FROM LEAF TOWARDS ROOT AND ADJUST MBR AND LHV VALUES*/
void adjustLHV(NODE parentNode)
{
    if (parentNode == NULL)
    {
        return;
    }
    for (int i = 0; i < parentNode->u.non_leaf_node.num_entries; i++)
    {
        parentNode->u.non_leaf_node.entries[i].lhv = calculateLHV(parentNode->u.non_leaf_node.entries[i]);
    }
    adjustLHV(parentNode->parent_ptr);
}
void adjustMBR(NODE parentNode)
{
    if (parentNode == NULL)
    {
        return;
    }
    for (int i = 0; i < parentNode->u.non_leaf_node.num_entries; i++)
    {
        parentNode->u.non_leaf_node.entries[i].mbr = calculateEntryMBR(parentNode->u.non_leaf_node.entries[i]);
    }
    adjustMBR(parentNode->parent_ptr);
}
void AdjustTree(NODE N, NODE NN, NODE *S, int s_size)
{
    // STOP IF ROOT LEVEL REACHED
    NODE Np = N->parent_ptr;
    NODE new_node = NULL;
    // PARENT = NULL; ROOT LEVEL
    if (!Np)
    {
        return;
    }
    // INSERT SPLIT NODE INTO PARENT
    if (NN)
    {
        // INSERT IN CORRECCT ORDER IF ROOM IN PARENT NODE
        if (Np->u.non_leaf_node.num_entries < MAX_CHILDREN)
        {
            InsertNode(Np, NN);
        }
        else
        {
            new_node = HandleOverFlowNode(Np, NN);
        }
        // Insert(Np, NN);
    }
    // ADJUST MBR AND LHV IN PARENT LEVEL
    // P = SET OF PARENT NODES FOR NODES IN S
    NODE *P = (NODE *)malloc(sizeof(struct Node) * MAX_POINTS);
    int numParents = 0;
    for (int i = 0; i < s_size; i++)
    {
        NODE parent = S[i]->parent_ptr;
        if (parent && !isInArray(P, numParents, parent))
        {

            P[numParents++] = parent;
        }
    }
    // ADJUST MBR AND LHV VALUES
    for (int i = 0; i < numParents; i++)
    {
        NODE parent = P[i];
        adjustMBR(parent);
        adjustLHV(parent);
    }

    // NEXT LEVEL
    AdjustTree(Np, new_node, P, numParents);
}

NODE ChooseLeaf(NODE n, Rectangle r, int h)
{
    /* RETURNS THE LEAF NODE IN WHICH TO PLACE A NEW RECTANGE*/

    /* IF N IS A LEAF, RETURN N*/
    if (n->is_leaf == 1)
    {
        return n;
    }

    /* IF N IS A NON-LEAF NODE, CHOOSE ENTRY (R, PTR, LHV) WITH MINIMUM LHV
    GREATER THAN H*/
    float min_LHV = INFINITY;
    NODE next_node = NULL;
    for (int i = 0; i < n->u.non_leaf_node.num_entries; i++)
    {
        if (n->u.non_leaf_node.entries[i].lhv > h && n->u.non_leaf_node.entries[i].lhv < min_LHV)
        {
            min_LHV = n->u.non_leaf_node.entries[i].lhv;
            next_node = n->u.non_leaf_node.entries[i].child_ptr;
        }
    }
    // IF ALL CHILDREN HAVE LHV LESS THAN H
    if (next_node == NULL)
    {
        // CHOOSE THE CHILD NODE WITH LARGEST LHV
        for (int i = 0; i < n->u.non_leaf_node.num_entries; i++)
        {

            if (n->u.non_leaf_node.entries[i].lhv > min_LHV)
            {
                min_LHV = n->u.non_leaf_node.entries[i].lhv;
                next_node = n->u.non_leaf_node.entries[i].child_ptr;
            }
        }
    }

    /* DESCEND UNTIL A LEAF NODE IS REACHED*/
    return ChooseLeaf(next_node, r, h);
}
NODE HandleOverFlow(NODE n, Rectangle rectangle)
{
    // E = SET OF ALL ENTRIES FROM N AND S-1 COOPERATING SIBLINGS
    NODE *S = (NODE *)malloc(sizeof(NODE) * MAX_POINTS);
    for (int i = 0; i < MAX_POINTS; i++)
    {
        S[i] = NULL;
    }
    S[0] = n;
    int numSiblings = 0;
    NODE parentNode = n->parent_ptr;
    // GO FROM CURRENT NODE TO ROOT FINDING SIBLINGS
    // WITH LESS THAN MAXIMUM POINTERS
    while (parentNode)
    {
        int index = -1;
        for (int i = 0; i < parentNode->u.non_leaf_node.num_entries; i++)
        {
            if (parentNode->u.non_leaf_node.entries[i].child_ptr == n)
            {
                index = i;
                break;
            }
        }
        if (index > 0 && parentNode->u.non_leaf_node.entries[index - 1].child_ptr->u.leaf_node.num_entries < M)
        {
            S[++numSiblings] = parentNode->u.non_leaf_node.entries[index - 1].child_ptr;
        }
        else if (index == 0 && parentNode->u.non_leaf_node.entries[index + 1].child_ptr->u.leaf_node.num_entries < M)
        {
            S[++numSiblings] = parentNode->u.non_leaf_node.entries[index + 1].child_ptr;
        }
        S[++numSiblings] = parentNode;
        parentNode = parentNode->parent_ptr;
    }

    int num_entries = 0;
    Rectangle E[MAX_POINTS];

    for (int i = 0; i < numSiblings; i++)
    {
        if (S[i]->is_leaf == 1)
        {
            for (int j = 0; j < S[i]->u.leaf_node.num_entries; j++)
            {
                E[num_entries++] = S[i]->u.leaf_node.entries[j].mbr;
            }
        }
        else
        {
            for (int j = 0; j < S[i]->u.non_leaf_node.num_entries; j++)
            {
                E[num_entries++] = S[i]->u.non_leaf_node.entries[j].mbr;
            }
        }
    }
    // ADD  R TO E
    E[num_entries++] = rectangle;

    // CHECK IF ANY NODE IS NOT FULL
    bool allFull = true;
    for (int i = 0; i < numSiblings; i++)
    {
        if (S[i]->is_leaf == 1)
        {
            if (S[i]->u.leaf_node.num_entries < M)
            {
                allFull = false;
            }
        }
        else
        {
            if (S[i]->u.non_leaf_node.num_entries < M)
            {
                allFull = false;
            }
        }
    }

    // IF ATLEAST ONE OF SIBLINGS IS NOT FULL
    if (!allFull)
    {

        // DISTRIBUTE E EVENLY AMONG THE S NODES ACCORDING TO THE HILBERT VALUE
    }
    // IF ALL COOPERATING SIBLINGS ARE FULL
    else
    {
        // CREATE A NEW NODE NN
        NODE NN = (NODE)malloc(sizeof(struct Node));
        // DISTRIBUTE E EVENLY AMONG THE S+1 NODES ACCORDING TO THE HILBERT VALUE

        return NN;
    }
    return NULL;
}
void Insert(NODE root, Rectangle rectangle)
{
    NODE leafNode = ChooseLeaf(root, rectangle, rectangle.h);
    NODE newLeafNode = NULL;
    if (leafNode->u.leaf_node.num_entries < M)
    {
        // LEAF NODE HAS SPACE
        int i = 0;
        while (i < leafNode->u.leaf_node.num_entries && leafNode->u.leaf_node.entries[i].mbr.h < rectangle.h)
        {
            i++;
        }

        for (int j = leafNode->u.leaf_node.num_entries; j > i; j--)
        {
            leafNode->u.leaf_node.entries[j] = leafNode->u.leaf_node.entries[j - 1];
        }
        // INSERT ACCORDING TO HILBERT ORDER AND RETURN
        LeafEntry entry;
        entry.mbr = rectangle;
        entry.obj_id = ++CURRENT_ID;
        leafNode->u.leaf_node.entries[i] = entry;
        leafNode->u.leaf_node.num_entries++;
        return;
    }
    else
    {
        // LEAF NODE IS FULL
        newLeafNode = HandleOverFlow(leafNode, rectangle);
        // RETURNS THE NEW LEAF IF SPLIT WAS INEVITABLE
    }

    // Propogate changes upward

    // FORM A SET S CONTAINING L: COOPERATING SIBLINGS AND NEW LEAF (IF ANY)
    NODE *S = (NODE *)malloc(sizeof(NODE) * MAX_POINTS);
    int numSiblings = -1;
    for (int i = 0; i < MAX_POINTS; i++)
    {
        S[i] = NULL;
    }
    S[++numSiblings] = leafNode;
    if (newLeafNode)
    {
        S[++numSiblings] = newLeafNode;
    }

    NODE parentNode = leafNode->parent_ptr;
    // GO FROM CURRENT NODE TO ROOT FINDING SIBLINGS
    // WITH LESS THAN MAXIMUM POINTERS
    while (parentNode)
    {
        int index = -1;
        for (int i = 0; i < parentNode->u.non_leaf_node.num_entries; i++)
        {
            if (parentNode->u.non_leaf_node.entries[i].child_ptr == S[0])
            {
                index = i;
                break;
            }
        }
        if (index > 0 && parentNode->u.non_leaf_node.entries[index - 1].child_ptr->u.leaf_node.num_entries < M)
        {
            S[++numSiblings] = parentNode->u.non_leaf_node.entries[index - 1].child_ptr;
        }
        else if (index == 0 && parentNode->u.non_leaf_node.entries[index + 1].child_ptr->u.leaf_node.num_entries < M)
        {
            S[++numSiblings] = parentNode->u.non_leaf_node.entries[index + 1].child_ptr;
        }
        S[0] = parentNode;
        parentNode = parentNode->parent_ptr;
    }
    AdjustTree(leafNode, newLeafNode, S, numSiblings + 1);

    // IF NODE SPLIT CAUSED ROOT TO SPLIT, CREATEA NEWROOT WITH CHILDREN
    // AS RESULTING NODES
    //  Check if the root split
    if (S[0]->parent_ptr == NULL)
    {
        NODE newRoot = (NODE)malloc(sizeof(struct Node));
        newRoot->is_leaf = 0;
        newRoot->u.non_leaf_node.num_entries = 2;
        newRoot->parent_ptr = NULL;

        newRoot->u.non_leaf_node.entries[0].child_ptr = S[0];
        // newRoot->u.non_leaf_node.entries[0].mbr = (S[0]->u.);

        newRoot->u.non_leaf_node.entries[1].child_ptr = S[1];
        // newRoot->u.non_leaf_node.entries[1].mbr = (S[1]->u.);

        S[0]->parent_ptr = newRoot;
        S[1]->parent_ptr = newRoot;

        root = newRoot;
    }
}

// SEARCH ALGORITHM
//-> NONLEAF - THOSE WITH MBR INTERSECTING THE QUERY WINDOW W
//-> LEAF - THOSE WITH MBR INTERSECTING THE QUERY WINDOW W
/*ALL RECTANGLES THAT OVERLAP A SEARCH RECTANGLE*/
bool intersects(Rectangle r1, Rectangle r2)
{
    return !(
        r1.top_right.x < r2.bottom_left.x || r2.top_right.x < r1.bottom_left.x ||
        r1.top_right.y < r2.bottom_left.y || r2.top_right.y < r1.bottom_left.y);
}
void search(NODE root, Rectangle rectangle, LeafEntry *results)
{
    if (root->is_leaf == 1)
    {
        for (int i = 0; i < root->u.leaf_node.num_entries; i++)
        {
            if (intersects(root->u.leaf_node.entries[i].mbr, rectangle))
            {
                results[num_results++] = root->u.leaf_node.entries[i];
            }
        }
    }
    else
    {
        for (int i = 0; i < root->u.non_leaf_node.num_entries; i++)
        {
            if (intersects(root->u.non_leaf_node.entries[i].mbr, rectangle))
            {
                search(root->u.non_leaf_node.entries[i].child_ptr, rectangle, results);
            }
        }
    }
}
// handle overflow
int numberOfSiblings(NODE n)
{
    NODE *S = (NODE *)malloc(sizeof(NODE) * MAX_POINTS);
    for (int i = 0; i < MAX_POINTS; i++)
    {
        S[i] = NULL;
    }
    S[0] = n;
    int numSiblings = 0;
    NODE parentNode = n->parent_ptr;
    // GO FROM CURRENT NODE TO ROOT FINDING SIBLINGS
    // WITH LESS THAN MAXIMUM POINTERS
    while (parentNode)
    {
        int index = -1;
        for (int i = 0; i < parentNode->u.non_leaf_node.num_entries; i++)
        {
            if (parentNode->u.non_leaf_node.entries[i].child_ptr == n)
            {
                index = i;
                break;
            }
        }
        if (index > 0 && parentNode->u.non_leaf_node.entries[index - 1].child_ptr->u.leaf_node.num_entries < M)
        {
            S[++numSiblings] = parentNode->u.non_leaf_node.entries[index - 1].child_ptr;
        }
        else if (index == 0 && parentNode->u.non_leaf_node.entries[index + 1].child_ptr->u.leaf_node.num_entries < M)
        {
            S[++numSiblings] = parentNode->u.non_leaf_node.entries[index + 1].child_ptr;
        }
        S[++numSiblings] = parentNode;
        parentNode = parentNode->parent_ptr;
    }
    return numSiblings;
}

NODE *cooperatingSiblings(NODE n)
{
    // taking s=2
    NODE *S = (NODE *)malloc(sizeof(NODE) * MAX_POINTS);
    for (int i = 0; i < MAX_POINTS; i++)
    {
        S[i] = NULL;
    }
    S[0] = n;
    int numSiblings = 0;
    NODE parentNode = n->parent_ptr;
    // GO FROM CURRENT NODE TO ROOT FINDING SIBLINGS
    // WITH LESS THAN MAXIMUM POINTERS
    int index = -1;
    for (int i = 0; i < parentNode->u.non_leaf_node.num_entries; i++)
    {
        if (parentNode->u.non_leaf_node.entries[i].child_ptr == n)
        {
            index = i;
            break;
        }
    }
    if (index > 0)
    {
        S[++numSiblings] = parentNode->u.non_leaf_node.entries[index - 1].child_ptr;
    }
    if (index < parentNode->u.non_leaf_node.num_entries - 1)
    {
        S[++numSiblings] = parentNode->u.non_leaf_node.entries[index + 1].child_ptr;
    }
    return S;
}
/*FIND THE LEAF NODE CONTAINING A RECTANGLE R*/
NODE findLeaf(NODE root, Rectangle rectangle)
{
    if (root->is_leaf == 1)
    {
        return root;
    }

    /* IF NON LEAF; FIND OVERLAPPING ENTRY*/
    int i;
    for (i = 0; i < root->u.non_leaf_node.num_entries; i++)
    {
        if (intersects(root->u.non_leaf_node.entries[i].mbr, rectangle))
        {
            break;
        }
    }
    return findLeaf(root->u.non_leaf_node.entries[i].child_ptr, rectangle);
}

int find_entry_index(NODE n, Rectangle rectangle)
{
    int index = -1;
    for (int i = 0; i < n->u.leaf_node.num_entries; i++)
    {
        if (n->u.leaf_node.entries[i].mbr.bottom_left.x == rectangle.bottom_left.x &&
            n->u.leaf_node.entries[i].mbr.bottom_left.y == rectangle.bottom_left.y &&
            n->u.leaf_node.entries[i].mbr.top_right.x == rectangle.top_right.x &&
            n->u.leaf_node.entries[i].mbr.top_right.y == rectangle.top_right.y)
        {
            index = i;
            break;
        }
    }
    return index;
}

/*DELETE(RECTANGLE R)*/
void delete(NODE root, Rectangle rectangle)
{
    // D1. FIND THE HOST LEAF
    // FIND THE HOST LEAF: N = LEAF CONTAINING ENTRY WITH RECTANGLE
    NODE n = findLeaf(root, rectangle);
    NODE parentNode = n->parent_ptr;
    NODE *S = (NODE *)malloc(sizeof(NODE) * MAX_POINTS);
    int numSiblings = 0;

    // D2. DELETE R: REMOVE R FROM NODE N
    // FIND INDEX OF ENTRY
    int index = find_entry_index(n, rectangle);

    if (index == -1)
    {
        printf("ENTRY NOT FOUND");
        return;
    }

    // ENTRY FOUND; INDEX>=0 DELETE ENTRY AT INDEX I
    for (int i = index; i < n->u.leaf_node.num_entries - 1; i++)
    {
        n->u.leaf_node.entries[i] = n->u.leaf_node.entries[i + 1];
    }
    n->u.leaf_node.num_entries--;

    // D3. IF NODE UNDERFLOWS: LESS THAN m ENTRIES
    if (n->u.leaf_node.num_entries < MIN_CHILDREN)
    {

        numSiblings = numberOfSiblings(n); // SIZE OF S
        S = cooperatingSiblings(n);
        int num_borrowed = 0;
        bool allUnderflow = false;
        int i = 0;

        for (int i = 1; i < numSiblings; i++)
        {
            if (S[i]->u.leaf_node.num_entries > MIN_CHILDREN)
            {
                allUnderflow = true;
            }
        }
        if (!allUnderflow)
        {
            while (n->u.leaf_node.num_entries < MIN_CHILDREN && i < numSiblings)
            {
                while (S[i] != NULL && S[i]->u.leaf_node.num_entries > MIN_CHILDREN && n->u.leaf_node.num_entries < MIN_CHILDREN)
                {
                    // LOGIC TO TRANSFER/BORROW AN ENTRY
                    n->u.leaf_node.entries[n->u.leaf_node.num_entries].mbr = S[i]->u.leaf_node.entries[0].mbr;
                    n->u.leaf_node.entries[n->u.leaf_node.num_entries].obj_id = S[i]->u.leaf_node.entries[0].obj_id;
                    n->u.leaf_node.num_entries++;
                    // SHIFT ALL ENTRIES ONE POSITION TO LEFT
                    for (int j = 0; j < S[i]->u.leaf_node.num_entries - 1; j++)
                    {
                        S[i]->u.leaf_node.entries[j] = S[i]->u.leaf_node.entries[j + 1];
                    }
                    S[i]->u.leaf_node.num_entries--;
                    num_borrowed++;
                    break;
                }
                i++;
            }
        }
        if (allUnderflow)
        {
            for (int i = 0; i < n->u.leaf_node.num_entries; i++)
            {
                // INSERT THE ENTRIES TO A COOPERATING SIBLING
                S[1]->u.leaf_node.entries[S[i]->u.leaf_node.num_entries] = n->u.leaf_node.entries[i];
                S[1]->u.leaf_node.num_entries++;
            }
            parentNode = n->parent_ptr;
            free(n);

            // DELETE THE ENTRY COMPLETELY
            for (int i = index; i < parentNode->u.non_leaf_node.num_entries - 1; i++)
            {
                parentNode->u.non_leaf_node.entries[i] = parentNode->u.non_leaf_node.entries[i + 1];
            }
            parentNode->u.non_leaf_node.num_entries--;
            // UPDATE LHV AND MBR OF NODE AND ENTRY
            adjustLHV(parentNode);
            adjustMBR(parentNode);
        }
        // BORROW SOME ENTRIES FROM S COOPERATING SIBLINGS Having more than minimum
        // IF ALL SIBLINGS READS TO UNDERFLOW; MERGE S+1 TO S NODES
        // ADJUST THE RESULTING NODES
    }
    adjustLHV(parentNode);
    adjustMBR(parentNode);
    // ADJUST MBR AND LHV IN PARENT LEVELS
    // FORM A SET S CONTAINING L AND COOPERATING SIBLINGS [IF UNDERFLOW HAD OCCURED]
    AdjustTree(n, NULL, S, numSiblings);
    free(S);
}

Rectangle calculateEntryMBR(NonLeafEntry entry)
{
    // FOR EACH NON LEAF ENTRY; CALCULATE MBR FROM CHILD  NODES
    // FIND -> LOWEST X, LOWEST Y  AND HIGHEST X, HIGHEST Y
    Rectangle mbr;
    int low_x = INFINITY;
    int low_y = INFINITY;
    int high_x = -INFINITY;
    int high_y = -INFINITY;

    NODE next_node = entry.child_ptr;
    // IF LEAF POINTER; FIND THE COORDINATES FROM ENTRIES
    if (next_node->is_leaf == 1)
    {
        for (int i = 0; i < next_node->u.leaf_node.num_entries; i++)
        {
            Rectangle obj_mbr = next_node->u.leaf_node.entries[i].mbr;
            low_x = (obj_mbr.bottom_left.x < low_x) ? obj_mbr.bottom_left.x : low_x;
            low_y = (obj_mbr.bottom_left.y < low_y) ? obj_mbr.bottom_left.y : low_y;
            high_x = (obj_mbr.top_right.x > high_x) ? obj_mbr.top_right.x : high_x;
            high_y = (obj_mbr.top_right.y > high_y) ? obj_mbr.top_right.y : high_y;
        }
    }
    else
    {
        // NON LEAF NODE:
        for (int i = 0; i < next_node->u.non_leaf_node.num_entries; i++)
        {
            Rectangle child_mbr = next_node->u.non_leaf_node.entries[i].mbr;
            low_x = (child_mbr.bottom_left.x < low_x) ? child_mbr.bottom_left.x : low_x;
            low_y = (child_mbr.bottom_left.y < low_y) ? child_mbr.bottom_left.y : low_y;
            high_x = (child_mbr.top_right.x > high_x) ? child_mbr.top_right.x : high_x;
            high_y = (child_mbr.top_right.y > high_y) ? child_mbr.top_right.y : high_y;
        }
    }
    mbr.bottom_left.x = low_x;
    mbr.bottom_left.y = low_y;
    mbr.top_right.x = high_x;
    mbr.top_right.y = high_y;
    mbr.h = HilbertValue(mbr.bottom_left, mbr.top_right);
    return mbr;
}

int calculateLHV(NonLeafEntry entry)
{
    int max_h = -INFINITY;
    NODE node = entry.child_ptr;
    // CALCULATE MAXIMUM H OF NODE
    if (node->is_leaf == 1)
    {
        for (int i = 0; i < node->u.leaf_node.num_entries; i++)
        {
            if (node->u.leaf_node.entries[i].mbr.h > max_h)
            {
                max_h = node->u.leaf_node.entries[i].mbr.h;
            }
        }
        return max_h;
    }
    else
    {
        // NON LEAF CHILD NODE
        for (int i = 0; i < node->u.non_leaf_node.num_entries; i++)
        {
            if (node->u.non_leaf_node.entries[i].mbr.h > max_h)
            {
                max_h = node->u.non_leaf_node.entries[i].mbr.h;
            }
        }
        return max_h;
    }
}
int main()
{
    FILE *fp;
    Point points[MAX_POINTS];
    int num_points = 0;

    // Open the file containing the data points
    fp = fopen("data.txt", "r");
    if (fp == NULL)
    {
        printf("Error opening file.\n");
        return 1;
    }

    // Read the data points from the file
    while (fscanf(fp, "%d %d\n", &points[num_points].x, &points[num_points].y) == 2)
    {
        num_points++;
        printf("POINT  = %d %d\n", points[num_points - 1].x, points[num_points - 1].y);
    }

    // Close the file
    fclose(fp);

    // LETS ASSUME THE RECTANGLES INITIALLY ARE THE DATA POINTS
    // AND THEIR CENTRES ARE THE POINTS TOO

    // CREATE RECTANGLES
    Rectangle rectangles[MAX_POINTS];
    for (int i = 0; i < num_points; i++)
    {
        rectangles[i].bottom_left = points[i];
        rectangles[i].top_right = points[i];
        rectangles[i].h = HilbertValue(points[i], points[i]);
    }

    HilbertRTree *Rtree = new_hilbertRTree();
    Rtree->root = new_node(1);
    for (int i = 0; i < MAX_POINTS; i++)
    {
        Insert(Rtree->root, rectangles[i]); // INSERT(NODE ROOT, RECTANGLE R)
    }

    return 0;
}
// void getCooperatingSiblings(NODE node, NODE* siblings, int s) {
//     NODE parent = node->parent_ptr;
//     int i;
//     int index = -1;
//     // Find the index of the given node in its parent's child pointers array
//     for (i = 0; i < parent->u.non_leaf_node.num_entries; i++) {
//         if (parent->u.non_leaf_node.entries[i].child_ptr == node) {
//             index = i;
//             break;
//         }
//     }
//     if (index == -1) {
//         // Node is not a child of its parent
//         return;
//     }

//     int numSiblings = 0;
//     // Check the s-1 siblings to the left of the given node
//     for (i = index-1; i >= 0 && numSiblings < s-1; i--) {
//         if (parent->u.non_leaf_node.entries[i].child_ptr->u.leaf_node.num_entries < M) {
//             siblings[numSiblings++] = parent->u.non_leaf_node.entries[i].child_ptr;
//         }
//     }

//     // Check the s-1 siblings to the right of the given node
//     for (i = index+1; i < parent->u.non_leaf_node.num_entries && numSiblings < s-1; i++) {
//         if (parent->u.non_leaf_node.entries[i].child_ptr->u.leaf_node.num_entries < M) {
//             siblings[numSiblings++] = parent->u.non_leaf_node.entries[i].child_ptr;
//         }
//     }
// }
