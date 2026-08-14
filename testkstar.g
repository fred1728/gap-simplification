e:=1;;
n:=0;;
k:=1;;
L:=ContractibleGcomplex("SL3Zs");;

stardup := List([1..k+1], i -> []);;
Stabe:=L!.stabilizer(n,e);;
stardup[1]:=[[e,1,Stabe]];;
dim:=1;;

for f in [1..L!.dimension(n+dim)] do
    bnd := List(L!.boundary(n+dim, f),x-> [AbsoluteValue(x[1]),x[2]]);;
    for cell in bnd do
        Li:= Filtered(stardup[dim],x -> x[1]=cell[1]);;
        g:=L!.elts[cell[2]];;
        for cell2 in Li do
            Stabf:=L!.stabilizer(n+dim,f);;
            Stabcell:=L!.stabilizer(n+dim-1,cell[1]);;
            K1:=Intersection(Stabf,Stabcell);;
            K2:=RightTransversal(Stabcell,K1);;
            for l in K2 do
                g1:= Inverse(l)*Inverse(g);;
                Add(stardup[dim+1],[f,g1,Stabf]);;
            od;;
        od;;
    od;;
od;;

for dim in [2..k] do
    for f in [1..L!.dimension(n+dim)] do
        bnd := List(L!.boundary(n+dim, f),x-> [AbsoluteValue(x[1]),x[2]]);;
        for cell in bnd do
            Li:= Filtered(stardup[dim],x -> x[1]=cell[1]);;
            g:=L!.elts[cell[2]];;
            for cell2 in Li do
                h:=cell2[2];;
                Stabf:=L!.stabilizer(n+dim,f);;
                Stabcell:=L!.stabilizer(n+dim-1,cell[1]);;
                K1:=Intersection(Stabf,Stabcell);;
                K2:=RightTransversal(Stabcell,K1);;
                for l in K2 do
                    g1:= h*Inverse(l)*Inverse(g);;
                    Add(stardup[dim+1],[f,g1,Stabf]);;
                od;;
            od;;
        od;;
    od;;
od;;


R:=[];;
star := List([1..k+1], i -> []);;
star[1]:=[[e,1,Stabe]];;
for dim in [2..k+1] do
    for cell in stardup[dim] do
        T:=true;;
        for cell2 in Filtered(star[dim],x->x[1]=cell[1]) do
            Stabcell:=L!.stabilizer(n+dim-1,cell[1]);;
            if Inverse(cell[2])*cell2[2] in Stabcell then
                T:=false;;
                Add(R,[cell,cell2]);;
                break;;
            fi;;
        od;;
        if T then
            Add(star[dim],cell);;
        fi;;
    od;;
od;;