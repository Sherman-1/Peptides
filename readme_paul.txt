# The 8th of Feb 2024
# pdbstyle-2.08-X
for i in $(seq 1 8)
do
	echo $i
	wget https://scop.berkeley.edu/downloads/pdbstyle/pdbstyle-2.08-${i}.tgz
	tar -xzf pdbstyle-2.08-${i}.tgz
	mv pdbstyle-2.08 pdbstyle-2.08-${i}
done
# pdbstyle-2.08_tree.txt
ls pdbstyle-2.08-*/*/ | awk '/^pdbstyle-2.08/ {dir=$0;gsub(/:/,"",dir)} ; /.ent$/ {print dir$0}' > pdbstyle-2.08_tree.txt
#
# astral-scopedom-seqres-gd-all-2.08-stable.fa
wget https://scop.berkeley.edu/downloads/scopeseq-2.08/astral-scopedom-seqres-gd-all-2.08-stable.fa
# class_counts.tsv
awk '/>/ {letter=gensub(/(.).*/,"\\1","g",$2) ; print letter}' astral-scopedom-seqres-gd-all-2.08-stable.fa | sort | uniq -c > class_counts.tsv
# class_g.txt
awk '/>/ {letter=gensub(/(.).*/,"\\1","g",$2) ; if(letter=="g"){print $0}}' astral-scopedom-seqres-gd-all-2.08-stable.fa | sed "s/>//" > class_g.txt 
# class_g.faa
faSomeRecords astral-scopedom-seqres-gd-all-2.08-stable.fa class_g.txt class_g.faa
# class_g_size.tsv
faSize -detailed class_g.faa > class_g_size.tsv
# class_g_20_to_70.txt
awk -F"\t" '$2 >= 20 && $2 <= 70 {print $1}' class_g_size.tsv > class_g_20_to_70.txt
# class_g_20_to_70.faa
faSomeRecords class_g.faa class_g_20_to_70.txt class_g_20_to_70.faa
# class_g_20_to_70_uniq.faa
# 50 sequences out of 2084 contain a "X" because there is two chains
python /store/EQUIPES/BIM/MEMBERS/paul.roginski/Eukaryotes/SCRIPTS/fasta_group_identical_sequences.py class_g_20_to_70.faa > class_g_20_to_70_uniq.faa
# class_g_20_to_70_uniq_sid.txt
grep ">" class_g_20_to_70_uniq.faa | sed "s/>//" | awk '{print $1}' > class_g_20_to_70_uniq_sid.txt
# class_g_20_to_70_uniq_sid_path.txt
# 46 sid not found, exactly those not starting with a 'd'
grep -f class_g_20_to_70_uniq_sid.txt pdbstyle-2.08_tree.txt > class_g_20_to_70_uniq_sid_path.txt
# class_g_20_to_70_uniq_pdb/
mkdir class_g_20_to_70_uniq_pdb/
cat class_g_20_to_70_uniq_sid_path.txt | while read file ; do cp $file class_g_20_to_70_uniq_pdb/ ; done
#
# classes_a_b_c_d_e.txt
for letter in a b c d e; do awk -v curr_letter=${letter} '/>/ {letter=gensub(/(.).*/,"\\1","g",$2) ; if(letter==curr_letter){print $0}}' astral-scopedom-seqres-gd-all-2.08-stable.fa | sed "s/>//" ; done > classes_a_b_c_d_e.txt
# classes_a_b_c_d_e.faa
faSomeRecords astral-scopedom-seqres-gd-all-2.08-stable.fa classes_a_b_c_d_e.txt classes_a_b_c_d_e.faa 
# classes_a_b_c_d_e_size.tsv
faSize -detailed classes_a_b_c_d_e.faa > classes_a_b_c_d_e_size.tsv
# classes_a_b_c_d_e_20_to_70.txt
awk -F"\t" '$2 >= 20 && $2 <= 70 {print $1}' classes_a_b_c_d_e_size.tsv > classes_a_b_c_d_e_20_to_70.txt
# classes_a_b_c_d_e_20_to_70.faa
faSomeRecords classes_a_b_c_d_e.faa classes_a_b_c_d_e_20_to_70.txt classes_a_b_c_d_e_20_to_70.faa
# classes_a_b_c_d_e_20_to_70_uniq.faa
# only 7 "X" out of more than 3000 sequences
python /store/EQUIPES/BIM/MEMBERS/paul.roginski/Eukaryotes/SCRIPTS/fasta_group_identical_sequences.py classes_a_b_c_d_e_20_to_70.faa > classes_a_b_c_d_e_20_to_70_uniq.faa
# classes_a_b_c_d_e_20_to_70_uniq_sid.txt
grep ">" classes_a_b_c_d_e_20_to_70_uniq.faa | sed "s/>//" | awk '{print $1}' > classes_a_b_c_d_e_20_to_70_uniq_sid.txt
# classes_a_b_c_d_e_20_to_70_uniq_sid_path.txt
# 6 sid not found, exactly those not starting with a 'd'
grep -f classes_a_b_c_d_e_20_to_70_uniq_sid.txt pdbstyle-2.08_tree.txt > classes_a_b_c_d_e_20_to_70_uniq_sid_path.txt
# classes_a_b_c_d_e_20_to_70_uniq_pdb/
mkdir classes_a_b_c_d_e_20_to_70_uniq_pdb/
cat classes_a_b_c_d_e_20_to_70_uniq_sid_path.txt | while read file ; do cp $file classes_a_b_c_d_e_20_to_70_uniq_pdb/ ; done
