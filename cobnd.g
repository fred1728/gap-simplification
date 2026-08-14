dim:=PositionProperty(bnd,IsEmpty)-2;
cobnd:=[];; #Cobnd[n+1] contains the info on n-cells.
for n in [0..dim] do
  #k:=2*(n+1)+1;#k:=1+2^(n+1);
  cobnd[n+1]:=List(bnd[n+1],i->[0]);
  for j in [1..Length(bnd[n+2])] do
    b:=bnd[n+2][j];
    k:=Length(b);
    for i in b{[2..k]} do
      cobnd[n+1][i][1]:=cobnd[n+1][i][1]+1;
      Add(cobnd[n+1][i],j);
    od;
  od;
od;
cobnd[dim+1]:=List(bnd[dim+1],a->[0]);