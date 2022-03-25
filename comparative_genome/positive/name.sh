for i in `cat a`
do
	sed -i "s/|.* /   /g" ${i}.phy
done
