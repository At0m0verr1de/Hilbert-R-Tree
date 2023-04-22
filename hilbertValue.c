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
