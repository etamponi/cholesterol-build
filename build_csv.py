"""
Build csv datasets using the json2csv.pl script, labels information from
labels.txt and the json files in the json/ directory.

json2csv.pl script expects 2 parameters: a json file to convert and the label to
apply to every entry in the json file.

The json files to convert are in the json/ directories. They are named as
follows::

    Q8NBP7_[aminoacid_from]_[position]_[aminoacid_to].json

Q8NBP7 is the id of the protein, `aminoacid_from` is the original aminoacid in
the chain, `position` indicates where in the chain the mutation occurs,
and `aminoacid_to` is the *final* aminoacid found in the chain after the
mutation.

labels.txt is formatted as follows::

    [position],[aminoacid_from],[aminoacid_to],[label]

The role of this script is twofold:

* map entries in the labels.txt script to json files and call json2csv.pl with
  the correct file-label pair.

* create a manual 10-fold in which the mutations in the training set do not
  appear in the test set.

It is important to note that some entries in the labels.txt file cannot be
converted to csv because the json file contains no data: the script takes this
into account by reporting the empty json file in a log (errors.log).
"""
from subprocess import check_output, CalledProcessError
import subprocess

__author__ = 'Emanuele Tamponi <emanuele.tamponi@gmail.com>'

HEADER = """\
Binding,SProtFT0,SProtFT1,SProtFT2,SProtFT3,SProtFT4,SProtFT5,\
SProtFT6,SProtFT7,SProtFT8,SProtFT9,SProtFT10,SProtFT11,SProtFT12,\
Interface,Relaccess,Impact,HBonds,SPhobic,CPhim,BCharge,SSGeom,\
Voids,MLargest1,MLargest2,MLargest3,MLargest4,MLargest5,MLargest6,\
MLargest7,MLargest8,MLargest9,MLargest10,NLargest1,NLargest2,\
NLargest3,NLargest4,NLargest5,NLargest6,NLargest7,NLargest8,\
NLargest9,NLargest10,Clash,Glycine,Proline,CisPro,dataset"""

PREFIX = "Q8NBP7"
COMMAND = "./json2csv.pl"
FOLDS = 10
DATASET_NAME = "cholesterol"
LIMIT = 100


def greedy_split(elements, k):
    # Like children choose teammates: strongest first.
    # The weakest team choose first.
    lengths = [len(e) for e in elements]
    splits = [[] for _ in xrange(k)]
    splits_length = [0] * k
    while sum(lengths) > 0:
        # Find maximum and its index
        l_max = max(lengths)
        i_max = lengths.index(l_max)
        # Find minimum length split (will get the next element)
        i_min = splits_length.index(min(splits_length))
        splits[i_min].append(elements[i_max])
        splits_length[i_min] += l_max
        lengths[i_max] = 0
    return splits


def main_equal_splits():
    mutations = {"GOF": [], "LOF": []}

    with open("labels.txt") as f, open("errors.log", "w") as err:
        for line in f:
            line = line.strip()
            position, aminoacid_from, aminoacid_to, label = line.split(",")
            json_file_name = "json/{}_{}_{}_{}.json".format(
                PREFIX, aminoacid_from, position, aminoacid_to
            )
            try:
                output = check_output([COMMAND, json_file_name, label],
                                      stderr=subprocess.STDOUT)
                output = output.split("\n")[:-1]  # Last one is empty
                if len(output) > LIMIT:
                    output = output[:LIMIT]
                mutations[label].append(output)
            except CalledProcessError as e:
                message = "{} {}: {}".format(label, json_file_name, e.output)
                err.write(message)
                print message,

    splits = {"GOF": greedy_split(mutations["GOF"], FOLDS),
              "LOF": greedy_split(mutations["LOF"], FOLDS)}
    splits = [splits["GOF"][i] + splits["LOF"][i] for i in range(FOLDS)]
    splits = [reduce(lambda a, b: a + b, split, []) for split in splits]

    fake_line = "0," * splits[0][0].count(",") + "---"
    with open("datasets/{}.csv".format(DATASET_NAME), "w") as complete:
        complete.write(HEADER + "\n")
        for i in xrange(FOLDS):
            fold_prefix = "{}_{:02d}".format(DATASET_NAME, i+1)
            with open("datasets/{}_0trn.csv".format(fold_prefix), "w") as lf, \
                    open("datasets/{}_1tst.csv".format(fold_prefix), "w") as tf:
                lf.write(HEADER + "\n")
                tf.write(HEADER + "\n")
                for j in xrange(FOLDS):
                    f = lf if j != i else tf
                    for line in splits[j]:
                        if i == j:
                            complete.write(line + "\n")
                        f.write(line + "\n")
                    if i == j:
                        complete.write(fake_line + "\n")


def main_one_test_per_mutation():
    mutations = []

    with open("labels.txt") as f, \
            open("labels_ok.txt", "w") as l, \
            open("errors.log", "w") as err:
        for line in f:
            line = line.strip()
            position, aminoacid_from, aminoacid_to, label = line.split(",")
            json_file_name = "json/{}_{}_{}_{}.json".format(
                PREFIX, aminoacid_from, position, aminoacid_to
            )
            try:
                output = check_output([COMMAND, json_file_name, label],
                                      stderr=subprocess.STDOUT)
                output = output.split("\n")[:-1]  # Last one is empty
                if len(output) > LIMIT:
                    output = output[:LIMIT]
                mutations.append(output)
                l.write(line + "\n")
            except CalledProcessError as e:
                message = "{} {}: {}".format(label, json_file_name, e.output)
                err.write(message)
                print message,

    fake_line = "0," * mutations[0][0].count(",") + "---"
    with open("datasets/{}_split.csv".format(DATASET_NAME), "w") as complete:
        complete.write(HEADER + "\n")
        for mutation in mutations:
            for line in mutation:
                complete.write(line + "\n")
            complete.write(fake_line + "\n")


if __name__ == "__main__":
    main_one_test_per_mutation()
