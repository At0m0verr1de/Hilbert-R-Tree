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
}
