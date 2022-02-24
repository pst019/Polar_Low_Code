#!/bin/csh

set array=(`seq 1 30`)
#set array={1..30}

echo $array

foreach a ($array)
echo $a

end
