git rev-parse --is-inside-work-tree >& gitCheck.tmp
if [ "$(<gitCheck.tmp)" == "true" ] 
then
	#echo "#define STAR_VERSION \"`git tag; git rev-parse HEAD`\"" | tr '\n' ' ' > VERSION
	echo "#define STAR_VERSION \"`git tag`\"" | tr '\n' ' ' > VERSION
fi
rm gitCheck.tmp
