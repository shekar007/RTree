#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))
#define MIN_ENTRIES 1
#define MAX_ENTRIES 4
#define NUM_POINTS 10000
#define filename "dataBig.txt"

typedef struct rtree_node rtree_node;
typedef struct rtree rtree;

// int *num_rects = 0;
/*
1. n1.entries initialise
*/
typedef struct pair
{
    float first;
    float second;
} pair;

typedef struct rectangle
{
    float x_min; // Minimum x-coordinate of the rectangle
    float x_max; // Maximum x-coordinate of the rectangle
    float y_min; // Minimum y-coordinate of the rectangle
    float y_max; // Maximum y-coordinate of the rectangle
    float x_centre;
    float y_centre;

} rectangle;

struct rtree_node
{
    rtree *tree;
    int is_leaf;                    // Flag indicating whether the node is a leaf node
    int num_entries;                // Number of entries in the node
    rectangle entries[MAX_ENTRIES]; // Array of rectangles stored in the node
    struct rtree_node *parent;      // Pointer to the parent node
    struct rtree_node **children;   // Pointer array to child nodes (if the node is not a leaf)
};

struct rtree
{
    rtree_node *root; // To keep track of root in r tree
};

typedef struct slice
{
    rectangle **sliceEle;
    int numRects;
} slice;

typedef struct search_result
{
    rectangle *matches;
    int num_matches;
    int max_limit;
} search_result;

// Function prototypes
void search(rtree_node *node, rectangle query);
void setRoot(rtree_node *node, rtree *root);
rectangle *createRect(float x_min, float y_min, float x_max, float y_max);
rtree_node *create_node(int isleaf);
void add_rectangle_to_node(rtree_node *node, rectangle rect);
int compareX(const void *a, const void *b);
void sort_rectangles_by_x_centre(rectangle *rects[], int num_rects);
int my_ceil(double num);
slice *putRectInSlices(rectangle **sortedXRectArr, int n);
rectangle **sort_rectangles_by_y_centre(rectangle **rects, int num_rects);
int compareY(const void *a, const void *b);
void sortRectInSlicesByY(slice *slices, int n);
rtree_node **packSortedSlicesIntoNodes(slice *slices, int noOfSlices);
rtree *sort_tile_recursive_internal(rtree_node **nodes, int num);
rectangle get_mbr(rtree_node *node);
rtree *create_r_tree(rectangle *rects[]);
search_result *search_r_tree(rtree *tree, rectangle query);
void search_node(rtree_node *node, rectangle query, search_result *result);
pair containsFully(rtree_node *node, rectangle x);
bool isLeafFull(rtree_node *node);
void adjust(rtree_node *node);
void setTreeRoot(rtree *root);
void insert(rtree_node *node, rectangle x);
float getArea(rectangle x);
void insertfunction(rtree *root, rectangle x);
int check_overlap(rectangle r1, rectangle r2);
void preorderTraverse(rtree_node *node);
void traverse(rtree root);

int main()
{
    FILE *file;
    file = fopen(filename, "r");
    rectangle *rects[NUM_POINTS];
    for (int i = 0; i < NUM_POINTS; i++)
    {
        float x, y;
        fscanf(file, "%f %f\n", &x, &y);
        rects[i] = createRect(x, y, x, y);
    }
    fclose(file);
    // Create an r tree from leaf nodes
    // insert a rectangle in the r tree

    rtree *my_tree = create_r_tree(rects);

    // inserting a point
    // rectangle *abc = createRect(1, 2, 1, 2);
    // insertfunction(my_tree, *abc);
    // Searching for a query rectangle

    printf("Searching for all points in query rectangle\n");
    rectangle query = *(createRect(804484, 747133, 804484, 747133));
    search_result *result = search_r_tree(my_tree, query);
    printf("Number of matches: %d\n", result->num_matches);
    for (int i = 0; i < result->num_matches; i++)
    {
        printf("Entry %d: x_min=%f, x_max=%f, y_min=%f, y_max=%f\n", i, result->matches[i].x_min, result->matches[i].x_max,
               result->matches[i].y_min, result->matches[i].y_max);
    }
    free(result);
    printf("Traversing R tree:\n");
    traverse(*my_tree);
    printf("Traversal Successful\n");
}

rectangle *createRect(float x_min, float y_min, float x_max, float y_max)
{
    /*
    Creates and returns a rectangle's pointer stored in the heap
    */
    rectangle *rect = calloc(1, sizeof(rectangle));
    rect->x_min = x_min;
    rect->y_min = y_min;
    rect->x_max = x_max;
    rect->y_max = y_max;
    rect->x_centre = (x_min + x_max) / 2.0;
    rect->y_centre = (y_min + y_max) / 2.0;
    return rect;
}

rtree_node *create_node(int isleaf)
{
    /*
    Creates and returns a rtree_node
    */
    rtree_node *node = (rtree_node *)malloc(sizeof(rtree_node));
    node->is_leaf = isleaf;
    node->num_entries = 0;
    node->children = calloc(MAX_ENTRIES, sizeof(rectangle *));
    for (int i = 0; i < MAX_ENTRIES; i++)
    {
        node->children[i] = NULL;
    }
    node->parent = NULL;
    return node;
}

void add_rectangle_to_node(rtree_node *node, rectangle rect)
{
    /*
    Adds a rectangle to a node. Reject the change if node
    is full
    */
    if (node->num_entries == MAX_ENTRIES)
    {
        // Reject
        printf("Overflow detected and prevented\n");
        return;
    }
    node->entries[node->num_entries] = rect;
    node->num_entries++;
    return;
}

int compareX(const void *a, const void *b)
{
    /*
    Compare function for qsort()
    - Compares the x co-ordinate for the
    centre of the rectangle
    */
    rectangle *rectA = *(rectangle **)a;
    rectangle *rectB = *(rectangle **)b;
    return (rectA->x_centre > rectB->x_centre) - (rectA->x_centre < rectB->x_centre);
}

void sort_rectangles_by_x_centre(rectangle *rects[], int num_rects)
{
    /*
    Sorts the rectangles by x co-ordinate
    of the centre
    */
    qsort(rects, num_rects, sizeof(rectangle *), compareX);
}

int my_ceil(double num)
{
    /*
    Returns an integer after ceiling the decimal
    part of the function
    */
    int int_part = (int)num;              // truncate decimal part
    double decimal_part = num - int_part; // get decimal part

    if (decimal_part > 0)
    {
        return int_part + 1; // round up to nearest integer
    }
    else
    {
        return int_part; // already an integer
    }
}

slice *putRectInSlices(rectangle **sortedXRectArr, int n)
{
    int S = ceil(sqrt(ceil((double)NUM_POINTS / MAX_ENTRIES)));
    slice *slices = calloc(S, sizeof(slice));
    /*
    If only one slice, deal with it seperately
    */
    if (S == 1)
    {
        slices[0].sliceEle = calloc(n, sizeof(rectangle *));
        slices[0].sliceEle = sortedXRectArr;
        slices[0].numRects = n;
        return slices;
    }
    /*
    Otherwise, insert the sorted rectangles from the sortedXRectArr and ensure
    that there are atleast m entries in each slice.
    */
    int x = 0;
    for (int i = 0; i < S; i++)
    {
        /*
            If we are dealing with last two slices, then ensure m rectangles in each slice.
        */
        if (i == S - 2)
        {
            int left = n - x;
            int secLast = S * MAX_ENTRIES;
            if (secLast > left)
                secLast = left;
            int last = left - S * MAX_ENTRIES;
            if (last < 0)
                last = 0;
            if (left - (S * MAX_ENTRIES) < MIN_ENTRIES)
            {
                while (last < MIN_ENTRIES)
                {
                    secLast--;
                    last++;
                }
            }
            slices[S - 2].sliceEle = calloc(secLast, sizeof(rectangle *));
            slices[S - 2].numRects = secLast;
            slices[S - 1].sliceEle = calloc(last, sizeof(rectangle *));
            slices[S - 1].numRects = last;
            /*
            Enter rectangles in the last two slices such that both have atleast m rectangles.
            */
            for (int k = 0; k < secLast; k++)
            {
                slices[S - 2].sliceEle[k] = sortedXRectArr[x++];
            }
            for (int k = 0; k < last; k++)
            {
                slices[S - 1].sliceEle[k] = sortedXRectArr[x++];
            }
            break;
        }
        /*
            If we are not dealing with the last two slices, then fill the slice fully
        */
        else
        {
            slices[i].sliceEle = calloc(S * MAX_ENTRIES, sizeof(rectangle *));
            slices[i].numRects = S * MAX_ENTRIES;
            for (int j = 0; j < S * MAX_ENTRIES; j++)
            {
                if (x == n)
                    break;
                slices[i].sliceEle[j] = sortedXRectArr[x++];
            }
        }
    }
    return slices;
}

int compareY(const void *a, const void *b)
{
    rectangle *rectA = *(rectangle **)a;
    rectangle *rectB = *(rectangle **)b;
    return (rectA->y_centre > rectB->y_centre) - (rectA->y_centre < rectB->y_centre);
}

rectangle **sort_rectangles_by_y_centre(rectangle **rects, int num_rects)
{
    qsort(rects, num_rects, sizeof(rectangle *), compareY);
    return rects;
}

/*
    Pick each slice one by one and sort the rectangles in that slice by y_centre
*/

void sortRectInSlicesByY(slice *slices, int n)
{
    for (int i = 0; i < n; i++)
    {
        ;
        slices[i].sliceEle = sort_rectangles_by_y_centre(slices[i].sliceEle, slices[i].numRects);
    }
}

/*
    Now, pack the sorted slices into nodes.
*/

rtree_node **packSortedSlicesIntoNodes(slice *slices, int noOfSlices)
{
    rtree_node **leafNodes = calloc(ceil((float)NUM_POINTS / MAX_ENTRIES), sizeof(rtree_node *));
    for (int i = 0; i < ceil((float)NUM_POINTS / MAX_ENTRIES); i++)
    {
        leafNodes[i] = create_node(1);
    }
    rectangle **rectArray = calloc(NUM_POINTS, sizeof(rectangle *));
    int x = 0;
    /*
        If there is only one slice, then only one node is made. This case
        is dealt seperately
    */

    if (noOfSlices == 1)
    {
        for (int j = 0; j < NUM_POINTS; j++)
        {
            rectArray[x++] = slices[0].sliceEle[j];
        }
        for (int j = 0; j < NUM_POINTS; j++)
        {
            add_rectangle_to_node(leafNodes[0], *rectArray[j]);
        }
        return leafNodes;
    }

    /*
        If more than one slice, then all slices before the last two are
        packed into nodes in the for loop.

    */

    for (int i = 0; i < noOfSlices - 2; i++)
    {
        for (int j = 0; j < MAX_ENTRIES * noOfSlices; j++)
        {
            rectArray[x++] = slices[i].sliceEle[j];
        }
    }

    /*
        Now, dealing with the last two slices such that each slice has atleast
        m rectangles.
    */

    int sizeSecLast = 0;
    int sizeLast = 0;

    int S = ceil(sqrt(ceil((float)NUM_POINTS / MAX_ENTRIES)));
    int maxPointsInASlice = S * MAX_ENTRIES;
    int remainingPoints = NUM_POINTS - ((noOfSlices - 2) * maxPointsInASlice);
    sizeSecLast = remainingPoints;

    if (remainingPoints < maxPointsInASlice + MIN_ENTRIES)
    {
        sizeLast = MIN_ENTRIES;
        sizeSecLast = sizeSecLast - MIN_ENTRIES;
    }
    else
    {
        sizeSecLast = maxPointsInASlice;
        sizeLast = remainingPoints - sizeSecLast;
    }

    /*
        Now, we know the no. of rectangles in last and second last slice
        and according to that, the retangles are packed.

    */

    for (int i = 0; i < sizeSecLast; i++)
    {
        rectangle *s = slices[noOfSlices - 2].sliceEle[0];
        rectArray[x++] = slices[noOfSlices - 2].sliceEle[i];
    }
    for (int i = 0; i < sizeLast; i++)
    {
        rectArray[x++] = slices[noOfSlices - 1].sliceEle[i];
    }
    int noOfLeaves = ceil((float)NUM_POINTS / MAX_ENTRIES);
    x = 0;
    for (int i = 0; i < noOfLeaves - 2; i++)
    {
        for (int j = 0; j < MAX_ENTRIES; j++)
        {
            if (x == NUM_POINTS)
                break;
            add_rectangle_to_node(leafNodes[i], *rectArray[x++]);
        }
    }

    /*
        Each leaf node must have atleast m entries. So, the number
        of rectangles for the last two nodes are calculated.
    */

    int leftRects = NUM_POINTS - x;
    int rectsInSecLastLeaf = MAX_ENTRIES;
    int rectsInLastLeaf = leftRects - MAX_ENTRIES;
    if (leftRects < MAX_ENTRIES)
    {
        rectsInSecLastLeaf = leftRects;
        rectsInLastLeaf = 0;
    }

    while (rectsInLastLeaf < MIN_ENTRIES)
    {
        rectsInSecLastLeaf--;
        rectsInLastLeaf++;
    }

    /*
        Now, the rectangles are added to the last two nodes as per the above calculation
    */

    for (int i = 0; i < rectsInSecLastLeaf; i++)
    {
        add_rectangle_to_node(leafNodes[noOfLeaves - 2], *rectArray[x++]);
    }

    for (int i = 0; i < rectsInLastLeaf; i++)
    {
        add_rectangle_to_node(leafNodes[noOfLeaves - 1], *rectArray[x++]);
    }

    for (int i = 0; i < noOfSlices; i++)
    {
        free(slices[i].sliceEle);
    }
    free(slices);
    return leafNodes;
}

rtree *sort_tile_recursive_internal(rtree_node **nodes, int num)
{
    /*
        This function recurses on nodes created till
        a root node is formed
    */
    if (num == 1)
    {
        // Base case
        // 1) Create an rtree struct
        // 2) return the struct
        rtree *tree = (rtree *)malloc(sizeof(rtree));
        tree->root = nodes[0];
        return tree;
    }

    // Recursive case
    // First, calculate number of new nodes
    int num_of_new_nodes;
    if (num % MAX_ENTRIES == 0)
    {
        num_of_new_nodes = (int)num / MAX_ENTRIES;
    }
    else
    {
        num_of_new_nodes = (int)(num / MAX_ENTRIES) + 1;
    }

    // Then, get MBR for each node and store it in a rectangle array
    rectangle rects[num];
    for (int i = 0; i < num; i++)
    {
        rects[i] = get_mbr(nodes[i]);
    }

    // Create nodes by grouping MBRs of child nodes
    rtree_node **new_nodes = (rtree_node **)malloc(sizeof(rtree_node *) * num_of_new_nodes);
    int rects_added = 0;
    for (int i = 0; i < num_of_new_nodes; i++)
    {
        new_nodes[i] = create_node(0);
        for (int j = 0; j < MAX_ENTRIES; j++)
        {
            if (rects_added == num)
            {
                break;
            }
            add_rectangle_to_node(new_nodes[i], rects[rects_added]);
            new_nodes[i]->children[j] = nodes[rects_added];
            nodes[rects_added]->parent = new_nodes[i];
            rects_added++;
        }
    }

    // Recurse sort_tile_recursive_internal till root node is created
    return sort_tile_recursive_internal(new_nodes, num_of_new_nodes);
}

rectangle get_mbr(rtree_node *node)
{
    /*
      This function is used to calculate MBR
      of a Node
    */
    rectangle mbr;
    int num_entries = node->num_entries;
    int x_min, y_min, x_max, y_max;
    x_min = node->entries[0].x_min;
    x_max = node->entries[0].x_max;
    y_min = node->entries[0].y_min;
    y_max = node->entries[0].y_max;
    for (int i = 1; i < num_entries; i++)
    {
        if (node->entries[i].x_max > x_max)
        {
            x_max = node->entries[i].x_max;
        }
        if (node->entries[i].x_min < x_min)
        {
            x_min = node->entries[i].x_min;
        }
        if (node->entries[i].y_max > y_max)
        {
            y_max = node->entries[i].y_max;
        }
        if (node->entries[i].y_min < y_min)
        {
            y_min = node->entries[i].y_min;
        }
    }
    mbr.x_max = x_max;
    mbr.x_min = x_min;
    mbr.y_max = y_max;
    mbr.y_min = y_min;
    return mbr;
}

rtree *create_r_tree(rectangle *rects[])
{
    /*
        This function manufactures an R Tree from an array of rectangles
        using the STR method of creating an R Tree
    */

    // Prints the array of rectangles
    for (int i = 0; i < NUM_POINTS; i++)
    {
        printf("Rectangle %d: (%.1f, %.1f) - (%.1f, %.1f)\n",
               i, rects[i]->x_min, rects[i]->y_min, rects[i]->x_max, rects[i]->y_max);
    }

    printf("Sorted by x:\n");
    // Sorts the rectangles using their x_centre
    sort_rectangles_by_x_centre(rects, NUM_POINTS);

    // Prints the sorted rectangles
    for (int i = 0; i < NUM_POINTS; i++)
    {
        printf("%f %f %f %f \n", rects[i]->x_min, rects[i]->y_min, rects[i]->x_max, rects[i]->y_max);
    }

    // Puts the sorted rectangles into slices
    slice *slices = putRectInSlices(rects, NUM_POINTS);

    // Prints the rectangles in each slice
    printf("Putting into slices:\n");
    for (int i = 0; i < ceil(sqrt((float)NUM_POINTS / MAX_ENTRIES)); i++)
    {
        printf("Slice %d:\n", i + 1);
        for (int j = 0; j < slices[i].numRects; j++)
        {
            printf("(%f %f %f %f) ", slices[i].sliceEle[j]->x_min, slices[i].sliceEle[j]->x_max,
                   slices[i].sliceEle[j]->y_min, slices[i].sliceEle[j]->y_max);
        }
        printf("\n");
    }
    printf("\n");

    // Sorts the rectangles in each slice by their y_centre co-ordinate
    sortRectInSlicesByY(slices, my_ceil(sqrt(my_ceil((float)NUM_POINTS / MAX_ENTRIES))));

    // Prints the sorted slices
    printf("Sorting rectangles in slices by Y:\n");
    for (int i = 0; i < ceil(sqrt((float)NUM_POINTS / MAX_ENTRIES)); i++)
    {
        printf("Slice %d:\n", i + 1);
        for (int j = 0; j < slices[i].numRects; j++)
        {
            printf("(%f %f %f %f) ", slices[i].sliceEle[j]->x_min, slices[i].sliceEle[j]->x_max,
                   slices[i].sliceEle[j]->y_min, slices[i].sliceEle[j]->y_max);
        }
        printf("\n");
    }
    printf("\n");

    // Packs the slices into nodes
    rtree_node **leaves = packSortedSlicesIntoNodes(slices, my_ceil(sqrt(my_ceil((float)NUM_POINTS / MAX_ENTRIES))));

    // Prints out the leaf nodes created
    printf("Packing sorted slices into nodes:\n");
    for (int i = 0; i < ceil((float)NUM_POINTS / MAX_ENTRIES); i++)
    {
        printf("Node %d:\n", i + 1);
        for (int j = 0; j < leaves[i]->num_entries; j++)
        {
            printf("(%f %f %f %f)", (leaves[i]->entries[j]).x_min, (leaves[i]->entries[j]).x_max,
                   (leaves[i]->entries[j]).y_min, (leaves[i]->entries[j]).y_max);
        }
        printf("\n");
    }
    int num_nodes = ceil((float)NUM_POINTS / MAX_ENTRIES);

    // Creates an R Tree by creating internal nodes
    // from leaf nodes till a root node is created
    rtree *tree = sort_tile_recursive_internal(leaves, num_nodes);
    return tree;
}

search_result *search_r_tree(rtree *tree, rectangle query)
{
    /*
        Initiales searching an R Tree and
        uses the search_node() to find matches
        Stores the result in 'result'.
    */
    search_result *result = (search_result *)malloc(sizeof(search_result));
    result->num_matches = 0;
    result->max_limit = 10;
    result->matches = (rectangle *)malloc(sizeof(rectangle) * 10);
    // Start search
    search_node(tree->root, query, result);
    return result;
}

void search_node(rtree_node *node, rectangle query, search_result *result)
{
    /*
        This function searches each node matching the query
        and ultimately adds all points inside the
        search query in 'result'
    */
    for (int i = 0; i < node->num_entries; i++)
    {
        // For every entry in a node,
        // checks for overlap
        if (check_overlap(node->entries[i], query))
        {
            if (node->is_leaf)
            {
                // Adds point to result if node is a leaf
                if (result->max_limit == result->num_matches)
                {
                    result->max_limit += 10;
                    result->matches = realloc(result->matches, (sizeof(rectangle) * (result->max_limit)));
                }
                result->matches[result->num_matches] = node->entries[i];
                result->num_matches++;
            }
            else
            {
                // Otherwise, recurses on the child of the node
                // with the matching query
                search_node(node->children[i], query, result);
            }
        }
    }
}

int check_overlap(rectangle r1, rectangle r2)
{
    /*
      Checks for overlap between 2 rectangles
    */
    if (r1.x_min > r2.x_max || r2.x_min > r1.x_max)
    {
        return 0;
    }
    if (r1.y_min > r2.y_max || r2.y_min > r1.y_max)
    {
        return 0;
    }
    return 1;
}

void preorderTraverse(rtree_node *node)
{
    /*
        Used to recursivly traverse through a node.
        First prints the node details,
        then recurses through its children, if any
    */
    // preorder = NULL
    if (node == NULL)
    {
        return;
    }
    if (node->is_leaf)
    {
        printf("Node is a LEAF \n");
        printf("Number of entries: %d\n", node->num_entries);
        for (int i = 0; i < node->num_entries; i++)
        {
            rectangle r = node->entries[i];
            printf("Entry %d: X=%f, Y=%f\n", i, r.x_centre, r.y_centre);
        }
    }
    if (!node->is_leaf)
    {
        // (top right point and bottom left point of the MBR being printed)
        printf("Node is a NOT LEAF \n");
        printf("printing data for node: \n");
        rectangle x = get_mbr(node);
        printf("top right point: X: %f, Y: %f\n", x.x_max, x.y_max);
        printf("bottom left point: X: %f, Y: %f\n", x.x_min, x.y_min);
        for (int i = 0; i < node->num_entries; i++)
        {
            printf("Child %d:\n", i);
            preorderTraverse(node->children[i]);
        }
    }
}

void traverse(rtree root)
{
    /*
        Starter function for pre-order traversal
    */
    preorderTraverse(root.root);
}

pair containsFully(rtree_node *node, rectangle x)
{
    rectangle parent = get_mbr(node);
    if (parent.x_min <= x.x_min && x.x_max <= parent.x_max && parent.y_min <= x.y_min && x.y_max <= parent.y_max)
    {
        float parent_area = (parent.y_max - parent.y_min) * (parent.x_max - parent.x_min);
        pair p = {0, parent_area};
        return p;
    }
    else
    {
        float ymax = max(x.y_max, parent.y_max);
        float ymin = min(x.y_min, parent.y_min);
        float xmax = max(x.x_max, parent.x_max);
        float xmin = min(x.x_min, parent.x_min);
        float area = (ymax - ymin) * (xmax - xmin);
        float parent_area = (parent.y_max - parent.y_min) * (parent.x_max - parent.x_min);
        pair p = {area - parent_area, parent_area};
        return p;
    }
}
bool isLeafFull(rtree_node *node)
{
    if (node->num_entries == MAX_ENTRIES)
    {
        return true;
    }
    else
    {
        return false;
    }
}

float getArea(rectangle x)
{
    return (x.x_max - x.x_min) * (x.y_max - x.y_min);
}

void adjust(rtree_node *node)
{
    if (node == node->tree->root)
    {
        rtree_node *newRoot = (rtree_node *)malloc(sizeof(rtree_node));
        newRoot->children = (rtree_node **)malloc(sizeof(rtree_node) * MAX_ENTRIES);
        newRoot->parent = (rtree_node *)malloc(sizeof(rtree_node));
        newRoot->tree = (rtree *)malloc(sizeof(rtree));

        // initialise newroot

        rtree_node *n1 = (rtree_node *)malloc(sizeof(rtree_node));
        rtree_node *n2 = (rtree_node *)malloc(sizeof(rtree_node));
        rtree_node *l1 = (rtree_node *)malloc(sizeof(rtree_node));
        n1->children = (rtree_node **)malloc(sizeof(rtree_node) * MAX_ENTRIES);
        n2->children = (rtree_node **)malloc(sizeof(rtree_node) * MAX_ENTRIES);
        n1->parent = (rtree_node *)malloc(sizeof(rtree_node));
        n2->parent = (rtree_node *)malloc(sizeof(rtree_node));
        l1->children = (rtree_node **)malloc(sizeof(rtree_node) * MAX_ENTRIES);
        if (node->children[0]->is_leaf == 1)
        {
            l1->is_leaf = 1;
        }
        l1->num_entries = 0;
        l1->parent = node;
        rtree_node *l2 = (rtree_node *)malloc(sizeof(rtree_node));
        l2->children = (rtree_node **)malloc(sizeof(rtree_node) * MAX_ENTRIES);
        if (node->children[0]->is_leaf == 1)
        {
            l2->is_leaf = 1;
        }
        l2->num_entries = 0;
        l2->parent = node;
        for (int i = 0; i < node->num_entries; i++)
        {
            rectangle t = get_mbr(node->children[i]);
            float ar = getArea(t);
            if (ar > getArea(get_mbr(l1)))
            {
                l1 = node->children[i];
                l2 = l1;
            }
            else if (ar > getArea(get_mbr(l2)))
            {
                l2 = node->children[i];
            }
        }
        n1->children[0] = l1;
        n1->num_entries = 1;
        n2->num_entries = 1;
        n1->parent = newRoot;
        n1->tree = node->tree;
        n2->children[0] = l2;
        n2->parent = newRoot;
        n2->tree = node->tree;
        n2->parent = node->parent;
        for (int i = 0; i < node->num_entries; i++)
        {
            rtree_node *f = node->children[i];
            if (f == l1 || f == l2)
            {
                continue;
            }
            else
            {
                float initn1 = getArea(get_mbr(n1));
                float initn2 = getArea(get_mbr(n2));
                n1->children[n1->num_entries] = f;
                n1->num_entries++;
                rectangle x1 = get_mbr(n1);
                n2->children[n2->num_entries] = f;
                n2->num_entries++;
                rectangle x2 = get_mbr(n2);
                if (getArea(x2) - initn2 > getArea(x1) - initn1)
                {
                    n2->num_entries--;
                    f->parent = n1;
                    // f is now child of n1
                }
                else
                {
                    // f is now child of n2

                    n1->num_entries--;
                    f->parent = n2;
                }
            }
        }
        newRoot->children = (rtree_node **)malloc(sizeof(rtree_node) * MAX_ENTRIES);
        newRoot->is_leaf = 0;
        newRoot->num_entries = 2;
        newRoot->parent = NULL;
        newRoot->tree = node->tree;
        newRoot->tree->root = newRoot;
        newRoot->children[0] = n1;
        newRoot->children[1] = n2;
        return;
    }

    rtree_node *n1 = (rtree_node *)malloc(sizeof(rtree_node));
    rtree_node *n2 = (rtree_node *)malloc(sizeof(rtree_node));
    n1->children = (rtree_node **)malloc(sizeof(rtree_node) * MAX_ENTRIES);
    n2->children = (rtree_node **)malloc(sizeof(rtree_node) * MAX_ENTRIES);
    n1->parent = (rtree_node *)malloc(sizeof(rtree_node));
    n2->parent = (rtree_node *)malloc(sizeof(rtree_node));

    rtree_node *l1 = (rtree_node *)malloc(sizeof(rtree_node));
    l1->children = (rtree_node **)malloc(sizeof(rtree_node) * MAX_ENTRIES);

    if (node->children[0]->is_leaf == 1)
    {
        l1->is_leaf = 1;
    }
    l1->num_entries = 0;
    l1->parent = node;
    rtree_node *l2 = (rtree_node *)malloc(sizeof(rtree_node));
    l2->children = (rtree_node **)malloc(sizeof(rtree_node) * MAX_ENTRIES);
    if (node->children[0]->is_leaf == 1)
    {
        l2->is_leaf = 1;
    }
    l2->num_entries = 0;
    l2->parent = node;
    int a1 = 0;
    int a2 = 0;
    for (int i = 0; i < node->num_entries; i++)
    {
        rectangle t = get_mbr(node->children[i]);
        float ar = getArea(t);
        if (ar > getArea(get_mbr(l1)))
        {
            l2 = l1;
            l1 = node->children[i];
            a2 = a1;
            a1 = i;
            // initialise l1 to be of 0 mbr area
        }
        else if (ar > getArea(get_mbr(l2)))
        {
            l2 = node->children[i];
            a2 = i;
        }
    }
    n1->children[0] = l1;
    n1->num_entries = 1;
    n1->is_leaf = 0;
    n2->is_leaf = 0;

    n2->num_entries = 1;
    n1->parent = node->parent;
    n2->children[0] = l2;
    n2->parent = node->parent;
    for (int i = 0; i < node->num_entries; i++)
    {
        rtree_node *f = node->children[i];
        if (i == a1 || i == a2) // change
        {
            continue;
        }
        else
        {
            float initn1 = getArea(get_mbr(n1));
            float initn2 = getArea(get_mbr(n2));
            n1->children[n1->num_entries] = f;
            n1->num_entries++;
            rectangle x1 = get_mbr(n1);
            n2->children[n2->num_entries] = f;
            n2->num_entries++;
            rectangle x2 = get_mbr(n2);
            if (getArea(x2) - initn2 > getArea(x1) - initn1)
            {
                n2->num_entries--;
                f->parent = n1;
                // f is now child of n1
            }
            else
            {
                // f is now child of n2

                n1->num_entries--;
                f->parent = n2;
            }
        }
    }

    rtree_node *par = node->parent;
    rtree_node **arr = par->children;
    for (int i = 0; i < par->num_entries; i++)
    {
        if (arr[i] == node)
        {
            arr[i] = n1;
            n1->parent = par;
            break;
        }
    }
    arr[par->num_entries] = n2;
    n2->parent = par;
    par->num_entries++;
    if (par->num_entries > MAX_ENTRIES)
    {
        adjust(par);
    }
}
void insert(rtree_node *node, rectangle x)
{
    if (!node->is_leaf)
    {
        rtree_node *minDiff = NULL;
        float parent_area = FLT_MAX;
        float minarea = FLT_MAX;

        for (int i = 0; i < node->num_entries; i++)
        {
            pair p = containsFully(node->children[i], x);
            float diff = p.first;
            if (minarea > diff)
            {
                minarea = diff;
                parent_area = p.second;
                minDiff = node->children[i];
            }
            else if (minarea == diff)
            {
                if (p.second < parent_area)
                {
                    minDiff = node->children[i];
                    parent_area = p.second;
                    minarea = diff;
                }
            }
        }
        insert(minDiff, x);
    }
    if (node->is_leaf)
    {
        if (isLeafFull(node))
        {
            rtree_node *n1 = (rtree_node *)malloc(sizeof(rtree_node));
            n1->children = (rtree_node **)malloc(sizeof(rtree_node) * MAX_ENTRIES);
            n1->is_leaf = 1;
            n1->parent = node->parent;
            n1->tree = node->tree;
            n1->num_entries = 1;
            n1->entries[0].x_max = x.x_max;
            n1->entries[0].x_centre = x.x_centre;
            n1->entries[0].x_min = x.x_min;
            n1->entries[0].y_centre = x.y_centre;
            n1->entries[0].y_max = x.y_max;
            n1->entries[0].y_min = x.y_min;
            rtree_node *par = node->parent;
            rectangle w = n1->entries[0];
            par->children[par->num_entries] = n1;
            par->num_entries++;
            rtree_node *q = par->children[par->num_entries - 1];
            rectangle v = q->entries[0];

            if (par->num_entries > MAX_ENTRIES)
            {
                adjust(node->parent);
            }
            // done insertion
        }

        else
        {
            // insertion if leaf node not full
            node->num_entries++;
            node->entries[node->num_entries - 1] = x;
        }
    }
}
void insertfunction(rtree *root, rectangle x)
{
    rtree_node *n = root->root;
    insert(n, x);
}
void setTree(rtree_node *node, rtree *root)
{
    // preorder = NULL
    if (node == NULL)
    {
        return;
    }
    if (node->is_leaf)
    {
        node->tree = root;
    }
    if (!node->is_leaf)
    {

        node->tree = root;
        for (int i = 0; i < node->num_entries; i++)
        {
            setTree(node->children[i], root);
        }
    }
}
void search(rtree_node *node, rectangle x)
{
    if (node == NULL)
        return;
    if (node->is_leaf == 1)
    {

        for (int i = 0; i < node->num_entries; i++)
        {
            rectangle r = node->entries[i];
            if (r.x_max == x.x_max && r.x_min == x.x_min && r.y_max == x.y_max && r.y_min == x.y_min)
            {
                printf("Match(es) found\n");
            }
        }
    }
    if (!node->is_leaf)
    {

        for (int i = 0; i < node->num_entries; i++)
        {
            search(node->children[i], x);
        }
    }
}

void setRoot(rtree_node *node, rtree *root)
{
    // preorder = NULL
    if (node == NULL)
    {
        return;
    }
    if (node->is_leaf)
    {
        node->tree->root = root->root;
    }
    if (!node->is_leaf)
    {

        node->tree->root = root->root;
        for (int i = 0; i < node->num_entries; i++)
        {
            setRoot(node->children[i], root);
        }
    }
}

void setTreeRoot(rtree *root)
{
    setTree(root->root, root);
    setRoot(root->root, root);
}
