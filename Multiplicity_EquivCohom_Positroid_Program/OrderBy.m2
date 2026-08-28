OrderBy = (original,ordering) -> (

n:=length(original);
ordered={};
for i from 0 to (n-1) do ordered=append(ordered, original#(ordering#i)); 
	
ordered
)

