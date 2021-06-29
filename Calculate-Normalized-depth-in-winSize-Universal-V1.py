########################################################
# nimarafati@gmail.com	                               #
# Please cite the script by referencing to github      #
# repository 					       #
########################################################
import argparse 
from time import sleep
import sys
import os

parser = argparse.ArgumentParser(description = "This script generates widnwos based depth of coverage using samtools depth output; Please make sure to use -a when running samtools depth to generate depth for all positions even if there is no depth. Use -Q 20 to select reads with mappinf quality above 20.")
parser.add_argument("-depth", "-d", help = "Depth file generated by samtools depth")
parser.add_argument("-window", "-w", help = "Non-overlapping widnow size (default 1000 bp)", type = int)
parser.add_argument("-chrome_size", "-fai", help = "Genome index file generated by samtools faidx or a tab delimited file with chromosomes size (Chr1	Size) wihtout header.")
parser.add_argument("-sample_list", "-s", help = "provide a sample list to name the columns accordingly. Please note that the order in the list should be similar to what you use in samtools depth command bam-list.")
args = parser.parse_args()

#def calc_window_depth (record):


chrome_dictionary = {}
chrome_size_file = open(args.chrome_size, 'r')
for chrome_size_line in chrome_size_file:
	chrome_size_line = chrome_size_line.split()
	key = chrome_size_line[0]
	value = chrome_size_line[1]
	chrome_dictionary[key] = value



if args.depth is None or args.chrome_size is None:
	print("Please provide depth and chromosome_size (fai) files")
	exit()

# Set default value to 1000 bp
if args.window is None:
	args.window = 1000
if args.sample_list:
	with open(args.sample_list) as f:
		 sample_line = [s.strip() for s in f]
	sample_line = [i for i in sample_line if i] #Remove empty elements ofa list 

depth_file = open(args.depth, 'r') 
open = 1
cntr = 1
prv_chrome = ""
for line in depth_file:
	line=line.strip()
	line_list = line.split('\t')
	chrome = line_list[0]

	if open  == 1:
	# Create empty array to save the average depth
	# Print header 	
		if args.sample_list:
			print('Chr      Start   End'),
			print('\t'.join(map(str,sample_line)))
		else:
			print('Chr      Start   End'),
			numbers = range(1,len(line[2:len(line)])+1)
			for i in range(1,len(line[2:len(line)])+1):
				print('\t sample'+str(i))
		sums = [0] * (int(len(line_list)) - 2) #Create empty list to save the average depths
	        pos = int(line_list[1]) - 1
		open = 0
		prv_chrome = chrome
	#Check if the query is on the same chromosome
	if  prv_chrome == chrome:
		#Sum as long as the query is within set window and count the counter
		if cntr < args.window :
			sums = [sum(x) for x in zip(sums, list(map(int, line_list[2:len(line_list)])))] #Sum as long as coordinate is within window size
			cntr = cntr + 1 
		#print if there is last query of the window and reset the variables
		elif cntr == args.window:
			sums = [sum(x) for x in zip(sums, list(map(int, line_list[2:len(line_list)])))]
			sums = [round(s / args.window, 2) for s in sums] #Average over the window size
			# Print per window average depth 
			print(str(chrome) + '\t' + str(pos) + '\t' + str(int(pos) + args.window)  + '\t' ),
			print('\t'.join(map(str,sums)))	
			
			#Reset variables
			sums = [0] * (int(len(line_list)) - 2)
			sums = [sum(x) for x in zip(sums, list(map(int, line_list[2:len(line_list)])))]
			cntr = 1
			pos = int(line_list[1]) # - 1 
			prv_pos = int(pos) #+ args.window
	# If it is new chromosome regardless of being the same size of the set window or not print and reset the variables
	elif prv_chrome != chrome and prv_chrome != '' :
		sums = [round(s / cntr, 2) for s in sums] #Average over the window size
		print(str(prv_chrome) + '\t' + str(prv_pos) + '\t' + str(chrome_dictionary[prv_chrome])  + '\t' ),
		print('\t'.join(map(str,sums)))
	
		#Reset varialbes
		sums = [0] * (int(len(line_list)) - 2)
		sums = [sum(x) for x in zip(sums, list(map(int, line_list[2:len(line_list)])))]
		cntr = 2
		prv_chrome = line_list[0]
		pos = int(line_list[1]) - 1

