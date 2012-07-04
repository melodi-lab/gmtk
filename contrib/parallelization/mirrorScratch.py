#!/usr/bin/python
"""
%InsertOptionParserUsage%

@author: 	Arthur Kantor
@contact: 	akantorREMOVE_THIS@uiuc.edu
@copyright:	Arthur Kantor 2008
@license: 	GPL version 3
@date: 		11/18/2008
@version: 	0.9
"""

import sys, re, struct, os, time, threading, signal;
from optparse import OptionParser, OptionGroup
from subprocess import *

#to parse the .torrent files
from sha import *
from BitTorrent.bencode import *

def usage():
 return """
mirrorScratch.py [--help] [options]  -m <pathToMirror>
Use bittorrent to distribute data to compute nodes.
The imagined use is to mirror some directory in the scratch space of the head node to the scratch spaces of the compute nodes.

Prerequisites:  
rtorrent and screen  is installed on all the compute nodes and the head node
btmaketorrent is installed on the starting node (typically the head node) 
A bittorrent tracker is running somewhere (ROCKS already has it running on the head node for compute node installation)

There are still a couple of efficiency bugs in the tracker and rtorrent (or perhaps I just don't know how to 
use them properly), but it will probably get sorted out in future versions.
1) The tracker reannounces peers once every 1800 seconds, even though I ask it to reannounce every 60 seconds
So, for the first 30 minutes, the seed node is the only one sharing anything.
2) In rtorrent, the transfer list (downloading chunks list) seems to get filled up with completed chunks, and eventually starts hitting memory limits and dropping peers with File chunk write error: Cannot allocate memory.
"""

def makeParser():
 parser = OptionParser()
 parser.usage=usage()
 parser.add_option("-u", "--trackerUrl", type="string", default="http://ifp-32:7625/announce",
                  help="bittorrent tracker URL. (default: %default)")
 parser.add_option("-b", "--bitTorrentScriptsPath", type="string", default="/opt/rocks/bittorrent",
                  help="The path to bittorrent scripts. (default: %default)")
 parser.add_option("-r", "--rTorrentPath", type="string", default="/cworkspace/ifp-32-1/hasegawa/programs/bin.x86_64/rtorrent",
                  help="The path to rtorrent binary. (default: %default)")
 parser.add_option("-t", "--torrent_dir", type="string", default="/cworkspace/ifp-32-1/hasegawa/programs/scripts/cluster/torrents",
                  help="Where to store torrent files. (default: %default)")
 parser.add_option("-s", "--save_in", type="string", default="/scratch/hasegawa",
                  help="The path where to store the files downloaded from a torrent. (default: %default)")
 parser.add_option("-m", "--pathToMirror", type="string", 
                  help="The path to directory or file to mirror. Required argument.")
 parser.add_option("--parse_dir_interval", type="int", default="60",
                  help="How often to check the torrent_dir for changes, in seconds. (default: %default sec)")
 parser.add_option("--display_interval", type="int", default="60",
                  help="How often the swarm participants print status, in seconds. (default: %default sec)")
 parser.add_option("--download_rate", type="int", default="9500",
                  help="Global maximim download rate limit, in kb. Currently there seems to be a bug in rtorrent, when the download rate is too fast (faster than the disk can write?), you get a 'File chunk write error: Cannot allocate memory.' and peers get dropped, and not picked up again until the next tracker announce.  This affects performance. (default: %default kb)")
 parser.add_option("--shutdown", action="store_true", default=False,
                  help="When in node mode, shutdown btlaunchmany.py (default: %default)")
 intUse = OptionGroup(parser,"OPTIONS FOR INTERNAL USE")
 intUse.add_option("-n", "--nodeMode", action="store_true", default=False,
                  help="Run script in node mode, i.e. ensure that exactly one copy of btlaunchmany.py is running (default: False)")
 parser.add_option_group(intUse)
 
 #command to include command line help in auto-generated documentation
 
 return parser
     
def main(argv):

 cmdLine=" ".join(argv[1:])
 opts, args = parser.parse_args()
 
 if opts.pathToMirror == None and opts.nodeMode == None :
	parser.error("need -m <pathToMirror>")
 
 
 
 mirror = RocksMirror(opts, cmdLine)
 if opts.nodeMode:
 	sys.exit(mirror.checkNode())
 else:
 	mirror.checkNodes()
 	if not opts.shutdown:
 		mirror.startBroadcast()
 		mirror.joinBroadcast()
  
class Mirror:
 """Uses bittorrent to efficiently mirror a directory across nodes in a cluster."""
 args =()
 cmdLine=""
 torrentFile=""
 nodes=[]
  
 def __init__(self, args, cmdLine):
	 self.args=args
	 self.cmdLine = cmdLine
	 
 def checkNodes(self):
	 """ Checks that the bittorrent is ready to talk on the compute nodes (e.g. btlaunchmany is running)."""
	 return NotImplemented
	
 def checkNode(self):
	"""Start/stop a node, as necessary.
	When starting, make sure that there is exactly one copy of btlaunchmany.py running on each compute node. 
	"""
	status=0
	#check if it's running already
	procs=self.getProcs()
	if self.args.shutdown:
		#shutdown requested
		if len(procs) == 0:
			print "### OK: already not running."
		else:
			for p in procs:
				words = str(p).split()
				pid = words[1]
				os.kill(int(pid), signal.SIGTERM)
			if len(self.getProcs()) == 0:
				print "### OK: stopped successfully."
			else:
				print "### ERROR: could not stop successfully."
				status=-1
	else:
    	#launch as needed 
		if len(procs) == 0:
		
			#cmd = "%s/btlaunchmany.py --max_upload_rate 0 --saveas_style 2 --parse_dir_interval %d --display_interval %d --save_in %s --torrent_dir %s  < /dev/null " % (
            #    self.args.bitTorrentScriptsPath,
            #    self.args.parse_dir_interval,
            #    self.args.display_interval,
    		#	self.args.save_in,
    		#	self.args.torrent_dir)
			cmd = "screen -md %s -d %s -O download_rate=%d -O schedule=watch_directory,%d,%d,load_start=%s/*.torrent -O schedule=untied_directory,%d,%d,stop_untied=" % (
                self.args.rTorrentPath,
    			self.args.save_in,
				self.args.download_rate,
    			self.args.parse_dir_interval,
    			self.args.parse_dir_interval,
    			self.args.torrent_dir,
    			self.args.parse_dir_interval,
    			self.args.parse_dir_interval,
				)
			print cmd
			bt = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
			#l = bt.stdout.readline()
    		
			#give it a second to start, and check if it's running - would be nicer to actually read the output and detach
			time.sleep(1)
			procs = self.getProcs()
			if len(procs) > 0:
				print "### OK: btlaunchmany started."
			else:
				print "### ERROR: could not start btlaunchmany." 			
		elif len(procs) == 1:
			print "### OK: btlaunchmany already running."
		else:
			print "### ERROR: more than one instance of btlaunchmany is running."
			status = -1	
 	
 	return status
 
 def getProcs(self):
 	"""Returns a list of lines, each line, is the output of the ps command for the btlaunchmany.py process"""
	#out=Popen("ps aux | grep 'btlaunchmany\.py .* --torrent_dir * %s'" % (self.args.torrent_dir), shell=True,  stdout=PIPE).communicate()[0]
	out=Popen("ps aux | grep 'rtorrent .*load_start=%s'" % (self.args.torrent_dir), shell=True,  stdout=PIPE).communicate()[0]
	procs=out.splitlines()
	procs = filter(lambda x : 'bin/sh' not in x and 'grep' not in x and 'SCREEN' not in x, procs) #don't count the bin/sh launcher on the head node
	#print out
	return procs
 
 def startBroadcast(self):
	 """Launches the mirroring between compute nodes. """
	 os.chdir(self.args.save_in)
	 self.torrentFile = "%s/%s.torrent" % (self.args.torrent_dir,self.args.pathToMirror.replace('/', '_').replace('.','_'))
	 print "creating torrent for '%s' to be saved in '%s'. Torrent saved in '%s'" % (self.args.pathToMirror, self.args.save_in, self.torrentFile)
	 cmd = "%s/btmaketorrent.py --piece_size_pow2 21 --target %s '%s' %s" % (
        self.args.bitTorrentScriptsPath,
        self.torrentFile,
		self.args.trackerUrl,
		self.args.pathToMirror)
	 retcode=0
	 retcode = call(cmd, shell=True)
	 if retcode != 0:
	 	raise("could not create torrent.")

 def joinBroadcast(self):
	 """Waits until all compute nodes have received the data. startBroadcast must have been called first."""
	 hash = self.getHash()
 	 print "Waiting for any peer to request torrent with hash %s from tracker..." % (hash) 

 	 while self.getTorrentStatus(hash) == None: 
 	    #print "### ERROR: Tracker has no nodes sharing the torrent with the above hash. Shut everything down and run the following by hand to see the error."
 	    #raise Exception("Tracker has no nodes sharing the torrent with the hash %s" % hash)
 	 	pass
 	 status=self.getTorrentStatus(hash)
 	 
	 print "polling the tracker until all nodes are seeding the torrent with hash %s..." % (hash) 
 	 
 	 print "complete\tdownloading"
 	 while status[0] != len(self.nodes):
 	 	status=self.getTorrentStatus(hash)
		print "%d\t%d" % (status[0],status[1])
		time.sleep(self.args.display_interval)

	 os.remove(self.torrentFile)
	 print "mirroring complete"
 def getHash(self):
    """Returns the torrent hash.  
        Adapted from btshowmetainfo.py.
    """ 
    metainfo_file = open(self.torrentFile, 'rb')
    metainfo = bdecode(metainfo_file.read())
    metainfo_file.close()
    announce = metainfo['announce']
    info = metainfo['info']
    #print info
    info_hash = sha(bencode(info))
    return info_hash.hexdigest()

 def getTorrentStatus(self, hash):
 	 """Returns a (	complete, downloading, downloaded) touple for the torrent with the given hash.
 	    Returns None if the tracker doesn't know such a torrent.
 	 """
 	 
	 trackerStatusUrl=self.args.trackerUrl.split('announce')[0] 	
 	 statusLines=Popen("wget %s -O - -q" % trackerStatusUrl, shell=True,  stdout=PIPE, stderr=PIPE).communicate()[0].splitlines()
 	 torrentStatus = filter(lambda x : hash in x, statusLines)
 	 if not torrentStatus:
 	 	return None
 	 parsedStatus=[]
 	 for v in torrentStatus[0].split('<code>'):
 	 	parsedStatus.append(v[0:v.index('<')])
 	 del parsedStatus[0:2];
 	 #print parsedStatus
 	 return map (lambda x: int(x), parsedStatus)
############################################################################################################### 
class RocksMirror(Mirror):
 """A Mirror class where the nodes are part of the ROCKS cluster."""

 
 def checkNodes(self):
	""" Checks that the bittorrent is ready to talk on the compute nodes (e.g. btlaunchmany is running).
		Raises Exception if at least one of the nodes errored.
	"""
	f=os.popen("rocks list host membership | tail -n +2 | sed -e 's/:.*//'")
	for n in f.readlines():
		self.nodes.append(n.strip())
	f.close()
	#print self.nodes
	#self.nodes=[self.nodes[0]]
	
	scriptPath = os.path.join(sys.path[0], sys.argv[0])
	
	threads=[]
	for n in self.nodes:
		cmd="ssh %s '%s --nodeMode %s < /dev/null'" %(n,scriptPath, self.cmdLine)
		t=CheckNodeThread(cmd,n)
		t.start()
		threads.append(t)
		
	groupStatus = 0
	badThread=""
	for t in threads:
		t.join()
		if t.status !=0:
			badThread=t.getName()
			groupStatus = t.status
	if groupStatus != 0:
		raise Exception("problem starting/stopping btlaunchMany in at least thread %s" % badThread)
	
	

class CheckNodeThread(threading.Thread):
 status=-2
 cmd=""
 
 def __init__(self,cmd,name):
 	threading.Thread.__init__(self,name=name)
 	self.cmd=cmd
 	
 def run(self):
	out=Popen(self.cmd, shell=True,  stdout=PIPE).communicate()[0]
	print "%s: %s" % (threading.currentThread().getName(), out.strip())
	if "### OK" in out:
		self.status = 0
	else:
		self.status = -1
        

#the parser is used for generating documentation, so create it always, and augment __doc__ with usage info  
#This messes up epydoc a little, but allows us to keep a single version of documentation for all purposes
parser = makeParser()
__doc__ = __doc__.replace("%InsertOptionParserUsage%\n", parser.format_help())

if __name__ == "__main__":
    main(sys.argv)

