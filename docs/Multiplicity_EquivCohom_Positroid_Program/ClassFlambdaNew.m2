ClassFlambdaNew = (f,lambda) -> (


--load "/Users/jflueg/Documents/PositroidCountProgram/ClassFlambdaNew.tex"
--assumes f is entered as a usual bounded affine permutation with indexing from 1 to n
--lastnum=take(f,{n-1,n-1});
--lastnum=lastnum#0-n;
--lastnum={lastnum};
--restnum=take(f,{0,n-2});
--f=join(lastnum,restnum);

n=length(f);
inc = toList (1..n);
k=(sum(f)-sum(inc))/n;
k=lift(k,ZZ);
lambcomp=inc-set lambda;

--for i from 0 to (n-1) do f=replace(i,f#i-1,f); --so now f has Macaulay ordering from 0 to n-1

fmod=apply(f, i -> i%n);
for i from 0 to (n-1) do (if fmod#i<0 then fmod#i=(fmod#i+n)); 

zeros=positions(fmod,i -> i==0);
for i from 0 to (length(zeros)-1) do fmod=replace(zeros#i,n,fmod);

ordering=SortOrder(fmod);
switchinv=inversePermutation ordering;

fordered=OrderBy(f,ordering);
jumps=positions(fordered, i -> i >n);

mu= n:0;
mu= toList mu;
numjumps=length(jumps)-1;
for i from 0 to numjumps do mu=replace(jumps#i,1,mu) ;
ordermuinv=join(positions(mu, j -> j==1),positions(mu, j -> j==0));
ordermu=inversePermutation ordermuinv;
omega=OrderBy(mu,ordermuinv);              --this is (1,1,...,1,0,...,0)

winv=ordermu_switchinv;
w=ordering_ordermuinv;

part1=take(w,{0,k-1});
wa=sort part1;
ord1=SortOrder(part1);

part2=take(w,{k,n-1});
wb=sort part2;
ord2=SortOrder(part2);
for i from 0 to (length(ord2)-1) do ord2=replace(i,ord2#i+k,ord2);

wgrass=join(wa,wb);
wgrassinv=inversePermutation wgrass;
for i from 0 to (n-1) do wgrass=replace(i,wgrass#i+1,wgrass);
for i from 0 to (n-1) do wgrassinv=replace(i,wgrassinv#i+1,wgrassinv);

grassperm= inversePermutation join(ord1, ord2);
grassinv= join(ord1, ord2);


-- here's what we have so far:  Our input is the bounded affine permutation f.
-- f= A^{-1} * Y * A * switchinv, where switch=ordering, A=ordermu, Y=t_\omega where \omega=(1,1,...,1,0,...,0) (k ones, n-k zeros), and the *'s are compositions of permutations 
-- next denoting w^{-1}=A*switchinv, we calculated that with grassperm we could make w into a Grassmannian permutation wgrass:
-- w=wgrass*grassperm, so A*switchinv=w^{-1}=grassperm^{-1}*wgrass^{-1}, and finally rearranging:
-- grassperm * A * switchinv = wgrass^{-1}
-- In sum, f = A^{-1} * (Y * grassperm^{-1}) * [grassperm * A * switchinv], the thing in [ ] is wgrass^{-1}
-- Now, we want to bring Y in (Y * grassperm^{-1}) to the other side to get (Y * grassperm^{-1})= X * Y, then we will get:
-- f = [A^{-1} * X] Y * [grassperm * A * switchinv] where the first [ ] is u, and the second [ ] is wgrass^{-1}
-- in the next set of lines, we calculate X, using that X = Y * grassperm^{-1} * Y^{-1} 

lattertwo=grassinv-n*omega; 
subs=OrderBy(omega,grassinv);
X=lattertwo+n*subs; 
vee=ordermuinv_X;
for i from 0 to (n-1) do vee=replace(i,vee#i+1,vee);




finalpol=0;

upermsa=permutations lambda;
upermsb=permutations lambcomp;

for i from 0 to (length(upermsa)-1) do (for j from 0 to (length(upermsb)-1) do(

uperm=join(upermsa#i,upermsb#j);

--this creates the polynomial we divide each piece of the sum by:
denompoly=1;
for p from 0 to (k-1) do (for q from (p+1) to (k-1) do (denompoly=denompoly*(z_(uperm#p)-z_(uperm#q))) );
for p from k to (n-1) do (for q from (p+1) to (n-1) do (denompoly=denompoly*(z_(uperm#p)-z_(uperm#q))) );

-- We're trying to calculate [x_v^w]|_u = [x_v]|_u * [x^w]|_u = [x_v]|_u * w_0\cdot([x_{w_0*w}]|_{w_0*u})

-- First: [x_v]|_u
xvee=DoubleSchubert(vee);
mapvars1={};
for p from 0 to (n-1) do (mapvars1=join(mapvars1,{z_(uperm#p)})); 
yvars=z_1..z_n; 
yvars= toList yvars;
mapvars1=join(mapvars1,yvars);
map1=map(FinalRing,SchubRing,mapvars1);
xvee=map1(xvee);

--Next; w_0\cdot([x_{w_0*w}]|_{w_0*u})
w0=0..(n-1);
w0=toList w0;
w0=rsort(w0);

ones=n:1;
ones=toList ones;

wshift=wgrass-ones;
w0w=w0_wshift;
w0w=w0w+ones;

upermshift=uperm-ones;
w0u=w0_upermshift;
w0u=w0u+ones;

mapvars2={};
for p from 0 to (n-1) do (mapvars2=join(mapvars2,{z_(w0u#p)})); 
mapvars2=join(mapvars2,yvars);
map2=map(FinalRing,SchubRing,mapvars2);


xw0w=DoubleSchubert(w0w);          --computes [x_{w_0*w}]
xw0wu=map2(xw0w);             --does the restriction [x_{w_0*w}]|_{w_0*u}

--finally, need w0\cdot to act on the left of [x_{w_0*w}]|_{w_0*u}:
mapvars3={};
for p from 0 to (n-1) do (mapvars3=join(mapvars3,{z_(n-p)}));
map3=map(FinalRing,FinalRing,mapvars3);
xw0wu=map3(xw0wu); 

polpiece=xvee*xw0wu/denompoly;

finalpol=finalpol+polpiece;
));

finalpol=lift(finalpol,FinalRing);

mapvars4={};
for p from 1 to n do (if isSubset({p},lambda) then (mapvars4=join(mapvars4,{h})) else mapvars4=join(mapvars4,{0}));
map4=map(FinalFinalRing,FinalRing,mapvars4);
finalfinalpol=map4(finalpol); 

finalpol = factor finalpol;
(finalpol, finalfinalpol)


)
