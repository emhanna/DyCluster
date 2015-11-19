# Modified by E.M. Hanna
#!/usr/bin/python
# (c) 2012 Jose Berengueres, Nazar Zaki and Dmitry Efimov. 2012.AUG.27
# UAE University, FIT, Bioinformatics Lab more on @bioAE Twitter

from optparse import OptionParser
import os

def main():
    print "Starting PEWCC ... please wait .."
    
    parser1 = OptionParser(description='predicts proteins from interactionfile.txt using MC or PE_EMH methods')

    parser1.add_option('-o','--overlap',  default=0.7,
                       help = "complexes that overlap more than this value will be merged, default is 0.7" )
    parser1.add_option('-s','--minhubsize',  default=41,
                       help = " min hub size for seed complexes, default is 41" )
    parser1.add_option('-r','--removeLowWeightedEdges', default=0.1,
                       help = "edges with less than this vlaue will be removed, default is 0.1" )
    parser1.add_option('-a','--addback', default=0.5,
                       help = "add back proteins whose conenctions to hub above this value, default is 0.5" )
    parser1.add_option('-m','--mode',  default='PE_EMH',
                       help = "choose PE_EMH or CD clsutering method, default is PE_EMH" )
    parser1.add_option('-f', '--fileinteractions', default='nofileselected.txt',
                       help = "protein interactions file" )

    options, args = parser1.parse_args()
    
    if options.mode == "PE_EMH":
        import PE_EMH
        PE_EMH.main(options)
    else:
        import CD2 as CD
        CD.main()

if __name__ == "__main__":
    main()
