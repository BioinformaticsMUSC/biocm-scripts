import csv

filtered_set = set()

with open('/Users/bryanwgranger/Documents/bioCM/deleon/analysis/fb/young/filtered_feature_bc_matrix/barcodes.tsv', 'r') as filtered_whitelist:
    csv_reader = csv.reader(filtered_whitelist, delimiter='\t')
    for cell_barcode in csv_reader:
        filtered_set.add(cell_barcode[0].strip('-1').strip())

with open('mapped_ADT_V3.txt', 'w') as outfile:
    with open('/Users/bryanwgranger/Documents/bioCM/deleon/3M-february-2018.txt') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')

        for row in csv_reader:

            if row[0].strip() in filtered_set:
                outfile.write('{}\n'.format(row[1].strip()))
