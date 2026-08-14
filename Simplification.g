SimplifiedFaceLattice:=function(Y)
local W , a, b, Schritt, OnceSimplifiedRegularCWComplex;


#######################################################
#######################################################
OnceSimplifiedRegularCWComplex:=function(W)
local Y, JoinCells, d, d1, c, d2, c2, empty, n, x, i, b,  cobnd, bnd,  dim , bool, orien,  pos ;

if IsBound(W!.orientation) then
Y:=RegularCWComplex(1*W!.boundaries,1*W!.orientation);
else
Y:=RegularCWComplex(1*W!.boundaries);

fi;

if IsBound(Y!.orientation) then bool:=true;
orien:=1*Y!.orientation;

else bool:=false; fi;

bnd:=Y!.boundaries;
cobnd:=Y!.coboundaries;



###################################################
###################################################
JoinCells:=function(d1,n)
                                 #The n-th cell in dimension d=d1-1 is removed
                                 #assuming it has a coboundary of size 2.
local V1,V2,V3,cob, d, dminus1, a, b, d2,d3, m, s, t, pos, poss, u, perm, cnt, F, Y2;

##
##CHECK IF REMOVAL SHOULD TAKE PLACE
if not (cobnd[d1][n][1] =2 and bnd[d1][n][1]>0) then return false; fi;
d2:=d1+1;
d3:=d2+1;
cob:=1*cobnd[d1][n];
if not SortedList(cobnd[d2][cob[2]])= SortedList(cobnd[d2][cob[3]]) then return false; fi;
V1:=BoundaryOfRegularCWCell(Y,d1-1,n);
V2:=BoundaryOfRegularCWCell(Y,d1,cobnd[d1][n][2]);
V3:=BoundaryOfRegularCWCell(Y,d1,cobnd[d1][n][3]);
if not 
 1+Size(V1)=Size(Intersection(V2,V3))  
then return false; fi;
##
##CHECK DONE

##
##REDUCE COBOUNDARIES OF BOUNDARIES OF nTH CELL
if d1>1 then
d:=d1-1;
for m in bnd[d1][n]{[2..Length(bnd[d1][n])]} do
    t:=1*cobnd[d][m]{[2..Length(cobnd[d][m])]};
    poss:=Position(t,n);
    Remove(t,poss);
    cobnd[d][m]:=Concatenation([Length(t)],t);
od;
fi;
##
##COBOUNDARIES OF BOUNDARIES REDUCED

##
##REMOVE nTH CELL, ITS COBOUNDARY, AND ADJUST ITS PRESENCE IN BOUNDARIES
##OF ITS COBOUNDARIES
bnd[d1][n]:=[0];
cobnd[d1][n]:=[0];

s:=bnd[d2][cob[2]];
s:=s{[2..Length(s)]};
pos:=Position(s,n);
Remove(s,pos);

if bool then a:=orien[d2][cob[2]][pos]; Remove(orien[d2][cob[2]],pos); fi;

t:=bnd[d2][cob[3]];
t:=t{[2..Length(t)]};
pos:=Position(t,n);
Remove(t,pos);

if bool then b:=orien[d2][cob[3]][pos]; Remove(orien[d2][cob[3]],pos); fi;

u:=SortedList(Concatenation(s,t));

if IsEmpty(u) then
    bnd[d2][cob[2]]:=[1,0];;
else
    bnd[d2][cob[2]]:=Concatenation([Length(u)],u);;
fi;;

#bnd[d2][cob[2]]:=Concatenation([Length(u)],u);

if bool then orien[d2][cob[2]]:=Concatenation(orien[d2][cob[2]],-a*b*orien[d2][cob[3]]); fi;
##
##nTH CELL AND ITS PRESENCE REMOVED


##
##FOR SECOND CELL OF DIMENSION n+1  REDUCE THE BOUNDARIES OF ITS 
##COBOUNDARIES

for m in cobnd[d2][cob[3]]{[2..Length(cobnd[d2][cob[3]])]} do
t:=bnd[d3][m]{[2..Length(bnd[d3][m])]};

pos:=Position(t,cob[3]);
Remove(t,pos);
bnd[d3][m]:=Concatenation([Length(t)],t);
if bool then Remove(orien[d3][m],pos); fi;
od;


##
##SECOND CELL  BOUNDARIES OF ITS COBOUNDARIES REDUCED

##
##REMOVE PRESENCE OF SECOND CELL IN COBOUNDARIES OF ITS BOUNDARIES
for m in bnd[d2][cob[3]]{[2..Length(bnd[d2][cob[3]])]} do
if cobnd[d1][m][1]>0 then

t:=cobnd[d1][m]{[2..Length(cobnd[d1][m])]};

pos:=Position(t,cob[3]);
t[pos]:=cob[2];
t:=SSortedList(t);
cobnd[d1][m]:=Concatenation([Length(t)],t);

fi;
od;


bnd[d2][cob[3]]:=[0];

cobnd[d2][cob[3]]:=[0];
#return true;
#end;


###################################################
###################################################



###################################################
######SIMPLIFICATION STARTS########################

# for d in [0..Dimension(Y)-1] do
# d1:=d+1;

# for n in [1..Length(bnd[d1])] do
# if cobnd[d1][n][1] =2 and bnd[d1][n][1]>0 then
# Print("d1=",d1,"\n");
# Print("n=",n,"\n");
# Print(bnd,"\n");
# JoinCells(d1,n);
# Print(bnd,"\n");   
# fi;
# od;

# od;

######SIMPLIFICATION DONE##########################
###################################################


perm:=[];  ##perm[d][i] will be the new position of
           ##the old i-th cell of dimension d;

for d in [1..Length(bnd)] do
perm[d]:=[];
cnt:=0;
for n in [1..Length(bnd[d])] do
if bnd[d][n][1]=0 then cnt:=cnt+1;
else
perm[d][n]:=n-cnt;
fi;
od;
od;

for d in [1..Length(bnd)] do
F:=1*Filtered([1..Length(bnd[d])],i->not bnd[d][i][1]=0);
bnd[d]:=1*((bnd[d]){F});
if bool and IsBound(orien[d]) then orien[d]:=orien[d]{F}; fi;

if d>1 then
dminus1:=d-1;
for x in bnd[d] do
    if x[2]>0 then
        for i in [2..Length(x)] do
            x[i]:=1*perm[dminus1][x[i]];
        od;
    fi;
od;
fi;
od;

#if bool then return RegularCWComplex(bnd, orien);
#if false then return RegularCWComplex(bnd, orien);

#if all lower-dimensional cells have been deleted, move all the remaining cells downwards
while Length(bnd) > 0 and IsEmpty(bnd[1]) do
    Remove(bnd, 1);
od;

Y2:=RegularCWComplex(bnd);
cobnd:=Y2!.coboundaries;

return true;
end;

#for c2 in bnd[1] do
#    c2[2]:=0;
#od;

###################################################
######SIMPLIFICATION STARTS########################

for d in [0..Dimension(Y)-1] do
d1:=d+1;

for n in [1..Length(bnd[d1])] do
if cobnd[d1][n][1] =2 and bnd[d1][n][1]>0 then
Print("d1=",d1,"\n");
Print("n=",n,"\n");
Print(bnd,"\n");
JoinCells(d1,n);
Print(bnd,"\n");   
fi;
od;

od;



return RegularCWComplex(bnd);
#fi;
end;
#######################################################
#######################################################

Schritt:=1;
W:=OnceSimplifiedRegularCWComplex(Y);
#Print("Schritt ",Schritt,W!.boundaries);

a:=Size(Y);
b:=Size(W);



while a>b do
W:=OnceSimplifiedRegularCWComplex(W);
a:=b;
b:=Size(W);
Schritt:=Schritt+1;
#Print("Schritt ",Schritt,W!.boundaries);
od;

return W;
end;