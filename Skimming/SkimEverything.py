############################################################################
#     This will skim all the samples in a samples.dat file                 #
############################################################################
import os, sys, optparse, subprocess


def mkdirIfNotExists(directory):
    if not os.path.exists(directory):
        os.mkdir(directory)



def makeSkim(nameOfSample, isMC, location, friendlocation, destiny):

    nameOfTheFile = location + "/" + nameOfSample + "/treeProducerSusyEdge/tree.root"
    nameOfTheFriend = (friendlocation+'/evVarFriend_'+ nameOfSample +'.root' if '/afs' in friendlocation else location+'/'+friendlocation+'/evVarFriend_'+ nameOfSample +'.root')

    nameOfTheDestinyFile = destiny + "/" + nameOfSample + "/treeProducerSusyEdge/tree.root"
    nameOfTheDestinyFriendFile = destiny + "/friends/evVarFriend_" + nameOfSample + ".root"
    mkdirIfNotExists(destiny)
    mkdirIfNotExists(destiny + "/" + nameOfSample)
    mkdirIfNotExists(destiny + "/" + nameOfSample + "/treeProducerSusyEdge/")
    mkdirIfNotExists(destiny + "/friends/")
    
    command = "./Skim " + nameOfTheFile + " " + nameOfTheFriend + " " + nameOfTheDestinyFile + " " + nameOfTheDestinyFriendFile + " " + isMC

    print command
    subprocess.call(command, shell=True)


if __name__ == "__main__":

    print 'Starting the skimming all the files...'
    parser = optparse.OptionParser(usage="usage: %prog [opts] -s FilenameWithSamples", version="%prog 1.0")
    #No options yet
    parser.add_option("-d", "--destiny", action="store", dest="destiny", default=".", help="Destiny directory")
    parser.add_option("-s", "--sample", action="store", dest="sample", default="./samples.dat", help="Sample file")

    (options, args) = parser.parse_args()
    originalSampleFile = options.sample
    destinySampleFile = "./Skimmed_samples.dat"
    destinyDirectory = options.destiny
    if os.path.exists(destinyDirectory):
        print "The destiny directory exits. It would be safer to check whether it is empty."
        sys.exit() 

    theDestinyFile = open(destinySampleFile, "w")




    for line in open(originalSampleFile).readlines():
        if(line[0] == "#"):
    	    theDestinyFile.write(line)
        else: 
            splitedLine = line.split()
            sourceOfFile = splitedLine[4]
            sourceOfFriendFile = splitedLine[5]
            splitedLine[4] = destinyDirectory + "/"
            splitedLine[5] = destinyDirectory + "/friends/"
            lineToWrite = ""
            BlankSpace = 0 
            for text in splitedLine:
              if BlankSpace == 0:
                lineToWrite = lineToWrite + text
                BlankSpace = 1
              else: 
                lineToWrite = lineToWrite + "    " + text
            lineToWrite = lineToWrite + "\n"
    	    theDestinyFile.write(lineToWrite)
            makeSkim(splitedLine[2], splitedLine[7], sourceOfFile, sourceOfFriendFile, destinyDirectory)


    theDestinyFile.close()
           
	





