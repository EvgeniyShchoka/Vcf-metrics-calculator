import json
import sys

def calculate_vcf_metrics(invcf, outvcf, outjson):
    """
    # Add AF (allele fraction, AD/DP) in the output.
    # Count variant per chromosome.
    # Calculate the ratio of heterozygous to homozygous variants. 
    # Calculate the transition to transversion ratio (Ti/Tv ratio).
    # Calculate the Allele Frequency (AF) for each variant and add it as a new field in the output file.
    # Remove variants with "0/0" genotypes from the output file. 
    # Create a new VCF file with AF fields. 
    # ##fileformat=VCFv4.2

    # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO           FORMAT  SAMPLE
    # chr1    100     rs1     A       T       30      PASS    AD=10,5;DP=15;AF=0.33  GT      0/1
    # chr1    150     rs2     C       G       45      PASS    AD=1,19;DP=20;AF=0.95  GT      0/1
    # chr1    150     rs2     C       A       45      PASS    AD=20,0;DP=20;AF=0  GT      0/0

    """

    with open (invcf, 'r') as file:
        # create empty dictionary to store a number of chr appearances
        chr_dict = {}
        # initialize variables to store the number of heterozygous and homozygous variants
        het = 0
        hom = 0
        # initialize variables to store the number of transitions and transversions
        transition = 0
        transversion = 0

        # read a conternt of the file string by string
        with open(outvcf, 'w') as out:
            for line in file:
                # chech whether this line is a heading
                if line.startswith("#"):
                    out.write(line)
                    continue
                
                # split line to separate strings  
                line_splitted = line.split("\t")

                # select INFO column
                info = line_splitted[7]

                # check whether this string contains '0/0' genotype
                genotype = line_splitted[9]
                if genotype == "0/0":
                    continue

                # split INfo to sseparate values and extract AD and AP
                info_list = info.split(";")
                
                AD = info_list[0].split(",")[1]
                DP = info_list[1].split("=")[1]

                # calculate AF
                AF = int(AD)/int(DP)

                # add value to the table
                AF_str = f";AF={AF}"
                line_splitted[7] = f"{info}{AF_str}"
                line_new = "\t".join(line_splitted)

                # write updated string to a new file
                out.write(line_new)



                # extract chr value
                chr = line_splitted[0]

                # checf if chr was previously else add it to dictionary
                if chr in chr_dict:
                    chr_dict[chr] += 1
                else:
                    chr_dict[chr] = 1



                # count a number of het and hom
                genotype = line_splitted[9]
                if "1/1" in genotype:
                    hom += 1
                elif "0/1" in genotype or "1/0" in genotype:
                    het += 1



                # extract REF and ALT
                ref = line_splitted[3]
                alt = line_splitted[4]

                transitions_ref = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

                if (ref, alt) in transitions_ref:
                    transition += 1
                else:
                    transversion += 1



        # count the ratio of het/hom
        het_hom_ratio = het / hom if hom != 0 else "N/A"

        # count the ratio of Ti/Tv
        ti_tv_ratio = transition/transversion if transversion != 0 else "N/A"



    data = {
        "number_of_variants": [
            {
                "chromosome": chr,
                "count": count
            } for chr, count in chr_dict.items()
        ],
        "het_hom_ratio": het_hom_ratio,  # number-of het / number-of hom
        "ti_tv_ratio": ti_tv_ratio,  # (A<>G or C<>T) / (A<>C, A<>T, G<>C, G<>T)
    }

    # save it to .json file
    with open(outjson, 'w') as out:
        json.dump(data, out, indent=4)



def main():
    # check if the number of arguments is correct
    if len(sys.argv) != 4:
        print("Usage: python3 generate_vcf_metric.py input.vcf output.vcf output.json")
        sys.exit(1)
    
    # store args in variables
    invcf = sys.argv[1]
    outvcf = sys.argv[2]
    outjson = sys.argv[3]

    # call main function
    calculate_vcf_metrics(invcf, outvcf, outjson)



if __name__ == "__main__":
    main()
