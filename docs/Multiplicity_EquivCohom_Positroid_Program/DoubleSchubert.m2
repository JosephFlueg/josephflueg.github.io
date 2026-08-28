DoubleSchubert = (perm) -> (

perm = toSequence perm;

if DS#?perm then DS#perm else (

len=length(perm);

if perm==(1..len) then output=1 else

(i:=len-1; 
while perm#i > perm#(i-1) do (i = i-1); 
i=i-1; 

ival=perm#i;

possj={};
for larg from i+1 to len-1 do (if perm#larg<ival then (possj=append(possj,perm#larg) ) ); 

j:=max(possj);

k:=i+1; 
while perm#i > perm#k and k<(len-1) do (k = k+1); 
if k==len-1 then (if perm#i < perm#k then k=k-1) else k=k-1;

newperm:=toList switch(i,k,perm);

permcover:={};
if i!=0 then for a from 0 to i-1 do(
if ((newperm#i>newperm#a) and isSubset((a..i),positions(newperm, m -> m<=newperm#a or m>=newperm#i)))
then permcover=append(permcover, switch(a,i,newperm)));

if length(permcover)>0 then (output=sum(permcover,DoubleSchubert)+(x_(i+1)-y_j)*DoubleSchubert(newperm)) else (output=(x_(i+1)-y_j)*DoubleSchubert(newperm))
);

DS#perm=output;
output
)
)
