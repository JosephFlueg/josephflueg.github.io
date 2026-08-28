SortOrder = (fmod) -> (

ordering:={};
temporary=fmod;

while length(temporary)>0 do(
minvalue= min temporary;
minpos=positions(fmod, i -> i==minvalue);
ordering=join(ordering,minpos);
temporary=select(temporary,i -> i!=minvalue)
);

ordering
)
