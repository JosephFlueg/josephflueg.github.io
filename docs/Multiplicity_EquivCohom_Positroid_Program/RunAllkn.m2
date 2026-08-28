RunAllkn = (k,n) -> (

k=lift(k,ZZ);
n=lift(n,ZZ);

smooth={};
singular={};
nonpoint={};

FinalFinalRing=QQ[h];
FinalRing=QQ[z_1..z_n];
SchubRing=QQ[x_1..x_n,y_1..y_n];

indc=1..n;
indc = toList indc;
onevec=k:1;
onevec=toList onevec;

ksubset=subsets(n,k);
allsitesw={};
uniqaffperm={};
uniqsitesw={};

DS= new MutableHashTable;

for ii from 0 to (length(ksubset)-1) do(
    increase=ksubset#ii;

    currperm=indc;
    --for jj from 0 to (k-1) do (currperm=replace(increase#jj,(currperm#jj)+n,currperm));
    currperm=apply(#currperm, jj->(if member(jj,set increase) then n else 0)) + currperm;

    combos=permutations currperm;

    for pp from 0 to (length(combos)-1) do (

	affperm=combos#pp;
	siteswap=affperm-indc;

	if min(siteswap)>=0 and max(siteswap)<=n then
	allsitesw=append(allsitesw,siteswap);
	);
);

 --this next section is Allen's code
scan(allsitesw, p->(nieuw := true; 
scan(n, i->(if member(rotate(i,p), uniqsitesw) then nieuw = false));
if nieuw then (uniqsitesw=append(uniqsitesw,p), uniqaffperm = append(uniqaffperm,p+indc));
));


for pp from 0 to (length(uniqsitesw)-1) do (
    siteswap=uniqsitesw#pp;
    affperm=siteswap+indc;
    
    for qq from 0 to (length(ksubset)-1) do(
	 
	lamb=ksubset#qq+onevec;
	(pol,hpol)=ClassFlambdaNew(affperm,lamb);

	hcoeff=sub(hpol,{h=>1});
	if hpol==1 or hcoeff==1 then smooth=append(smooth,{siteswap,lamb});
	if hpol==0 or hcoeff==0 then nonpoint=append(nonpoint,{siteswap,lamb});
	if hcoeff>1 then singular=append(singular,{siteswap,lamb,hcoeff});


)
) 
)



