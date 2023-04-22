int num_results = 0;
/*ALL RECTANGLES THAT OVERLAP A SEARCH RECTANGLE*/
void search(NODE root, Rectangle rectangle, LeafEntry* results){
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
