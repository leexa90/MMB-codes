for i in `ls -d TPP_FL_sln.out.*.pdb`; 
  do python make_last110pdb.py 2gis_A.pdb A $i C;
  ./MMB.2_17.Linux -c command121;
  cp trajectory.121.pdb ${i}trajectory121 ;
  cp last.121.pdb top${i}last121;
  python get_long_add_first_three.py top${i}last121 C;
  #rm top${i}last121;
  done


