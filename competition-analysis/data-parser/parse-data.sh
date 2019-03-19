echo "Parsing data files for static scenarios"
counter=0;
#for i in `find ../submissions/ -maxdepth 1 -mindepth 1 -type d -printf '%f\n' |sort`
for i in `find ../submissions/ -maxdepth 1 -mindepth 1 -type d | cut -f 4 -d'/' |sort`
do
	echo "Parsing: ${i}...";
	if [ ${counter} -eq 0 ];
	then
		./parse-static -r 50 -p 20 -D alldata-static.csv -P ../submissions/${i}/ -a ${i} -H
	else
		./parse-static -r 50 -p 20 -D alldata-static.csv -P ../submissions/${i}/ -a ${i}
	fi
	counter=$[${counter}+1];
done
echo "Parsing data files for dynamic scenarios"
counter=0;
#for i in `find ../submissions/ -maxdepth 1 -mindepth 1 -type d -printf '%f\n' |sort`
for i in `find ../submissions/ -maxdepth 1 -mindepth 1 -type d | cut -f 4 -d'/' |sort`
do
	echo "Parsing: ${i}...";
	if [ ${counter} -eq 0 ];
	then
		./parse-dynamic -r 50 -p 20 -D alldata-dynamic.csv -P ../submissions/${i}/ -a ${i} -H
	else
		./parse-dynamic -r 50 -p 20 -D alldata-dynamic.csv -P ../submissions/${i}/ -a ${i}
	fi
	counter=$[${counter}+1];
done
