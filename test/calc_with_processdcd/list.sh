awk '{i++;if($3=="CGA"){print i > "con_p1.list"};if($3=="CGB"){print i > "con_p2.list"};if($3=="CGC"){print i > "con_p3.list"}}' ../sys.pdb
