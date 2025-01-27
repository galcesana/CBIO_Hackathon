import FromVCFtoPresenceMatrix
import FromPresenceMatrixtoDistanceMatrix
import FromDistanceMatrixtoNewickFormat


if __name__ == "__main__":
    persons = ["PD4781","PD5117","PD5147","PD5163","PD5179","PD5182","PD5847","PD6629","PD6634","PD6646",
               "PD9478"]
    for person in persons:
        presence_path = "Results/"+person+"_presence.tsv"
        distance_path = "Results/"+person+"_distance.tsv"
        FromVCFtoPresenceMatrix.main("vcf/"+ person +".vcf", presence_path)
        FromPresenceMatrixtoDistanceMatrix.main(presence_path, distance_path)
        FromDistanceMatrixtoNewickFormat.main(distance_path, person+"_tree.nwk")

