if [ $# -lt 2 ] ; then
	echo "Usage $(basename $0) seed_lower_bound seed_upper_bound"
	exit
fi

min=$1
max=$2

cutoff=30000

i=$min

while [ $i -le $max ]
do
	./myMinWVC ../myMinWDS/bio-yeast $i $cutoff
	./myMinWVC ../myMinWDS/ca-AstroPh $i $cutoff
	./myMinWVC ../myMinWDS/ca-citeseer $i $cutoff
	./myMinWVC ../myMinWDS/ca-CSphd $i $cutoff

	i=$[$i+1]
done
