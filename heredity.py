import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):
                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.
    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    probability = 1
    #zero gene people and have not trait people
    zero_gene = set()
    not_have_trait = set()
    for i in people.keys():
        if i not in one_gene and i not in two_genes:
            zero_gene.add(i)
        if i not in have_trait:
            not_have_trait.add(i)

    #zero_trait_people
    for i in zero_gene:
        #if mother father data not there
        if people[i]["mother"] == None:
            prob = PROBS["gene"][0]
        #if mother father data found
        elif people[i]["mother"] != None:
            prob = float(1)
            mother = people[i]["mother"]
            father = people[i]["father"]

            if mother in zero_gene and father in zero_gene:
                prob *= (1-PROBS["mutation"]) * (1-PROBS["mutation"])
            elif mother in zero_gene and father in one_gene:
                prob *= (1-PROBS["mutation"]) * 0.5
            elif mother in zero_gene and father in two_genes:
                prob *= (1-PROBS["mutation"]) * PROBS["mutation"]

            elif mother in one_gene and father in zero_gene:
                prob *= 0.5 * (1-PROBS["mutation"])
            elif mother in one_gene and father in one_gene:
                prob *= 0.5 * 0.5
            elif mother in one_gene and father in two_genes:
                prob *= 0.5 * PROBS["mutation"]


            elif mother in two_genes and father in zero_gene:
                prob *= PROBS["mutation"] * (1-PROBS["mutation"])
            elif mother in two_genes and father in one_gene:
                prob *= PROBS["mutation"] * 0.5
            elif mother in two_genes and father in two_genes:
                prob *= PROBS["mutation"] * PROBS["mutation"]


        #if probability of have_trait
        if i in have_trait:
            prob *= PROBS["trait"][0][True]
        elif i in not_have_trait:
            prob *= PROBS["trait"][0][False]


        probability *= prob
    for i in one_gene:

        if people[i]["mother"] == None:
            prob = PROBS["gene"][1]

        elif people[i]["mother"] != None:
            prob = 1
            mother = people[i]["mother"]
            father = people[i]["father"]

            if mother in zero_gene and father in zero_gene:
                prob *= PROBS["mutation"] * (1-PROBS["mutation"]) + (1-PROBS["mutation"]) * PROBS["mutation"]
            elif mother in zero_gene and father in one_gene:
                prob *= PROBS["mutation"] * 0.5 + (1-PROBS["mutation"]) * 0.5
            elif mother in zero_gene and father in two_genes:
                prob *= PROBS["mutation"] * PROBS["mutation"] + (1-PROBS["mutation"]) * (1-PROBS["mutation"])

            elif mother in one_gene and father in zero_gene:
                prob *= 0.5 * (1-PROBS["mutation"]) + 0.5 * PROBS["mutation"]
            elif mother in one_gene and father in one_gene:
                prob *= 0.5 * 0.5 + 0.5 * 0.5
            elif mother in one_gene and father in two_genes:
                prob *= 0.5 * PROBS["mutation"] + 0.5 * (1-PROBS["mutation"])


            elif mother in two_genes and father in zero_gene:
                prob *= (1-PROBS["mutation"]) * (1-PROBS["mutation"]) + PROBS["mutation"] * PROBS["mutation"]
            elif mother in two_genes and father in one_gene:
                prob *= (1-PROBS["mutation"]) * 0.5 + PROBS["mutation"] * 0.5
            elif mother in two_genes and father in two_genes:
                prob *= (1-PROBS["mutation"]) * PROBS["mutation"] + PROBS["mutation"] * (1-PROBS["mutation"])
        if i in have_trait:
            prob *= PROBS["trait"][1][True]
        elif i in not_have_trait:
            prob *= PROBS["trait"][1][False]
        probability *= prob

    for i in two_genes:

        if people[i]["mother"] == None:
            prob = PROBS["gene"][2]

        elif people[i]["mother"] != None:
            prob = 1
            mother = people[i]["mother"]
            father = people[i]["father"]

            if mother in zero_gene and father in zero_gene:
                prob *= PROBS["mutation"] * PROBS["mutation"]
            elif mother in zero_gene and father in one_gene:
                prob *= PROBS["mutation"] * 0.5
            elif mother in zero_gene and father in two_genes:
                prob *= PROBS["mutation"] * (1-PROBS["mutation"])

            elif mother in one_gene and father in zero_gene:
                prob *= 0.5 * PROBS["mutation"]
            elif mother in one_gene and father in one_gene:
                prob *= 0.5 * 0.5
            elif mother in one_gene and father in two_genes:
                prob *= 0.5 * (1-PROBS["mutation"])


            elif mother in two_genes and father in zero_gene:
                prob *= (1-PROBS["mutation"]) * PROBS["mutation"]
            elif mother in two_genes and father in one_gene:
                prob *= (1-PROBS["mutation"]) * 0.5
            elif mother in two_genes and father in two_genes:
                prob *= (1-PROBS["mutation"]) * (1-PROBS["mutation"])
        if i in have_trait:
            prob *= PROBS["trait"][2][True]
        elif i in not_have_trait:
            prob *= PROBS["trait"][2][False]
        probability *= prob


    return probability


def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    for i in probabilities.keys():
        if i not in one_gene and i not in two_genes:
            gene = 0
        elif i in one_gene:
            gene = 1
        else:
            gene = 2
        probabilities[i]["gene"][gene] +=  p

        if i in have_trait:
            probabilities[i]["trait"][True] += p
        if i not in have_trait:
            probabilities[i]["trait"][False] += p



def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """

    for i in probabilities.keys():
        normal_gene = sum(probabilities[i]["gene"].values())
        normal_trait = sum(probabilities[i]["trait"].values())

        for j in probabilities[i]["gene"].keys():
            probabilities[i]["gene"][j] /= normal_gene

        for j in probabilities[i]["trait"].keys():
            probabilities[i]["trait"][j] /= normal_trait



if __name__ == "__main__":
    main()