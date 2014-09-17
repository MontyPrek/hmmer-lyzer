import csv
import sys
import matplotlib.pyplot as plt
import copy

class HmmerAnalyze:

    def __init__(s, sns_filename, hxt_filename, 
                 delimiter = ',', quotechar = '"'):
        s.sns_filename = sns_filename
        s.hxt_filename = hxt_filename

        s.delimiter = delimiter
        s.quotechar = quotechar
        
        s.sns_score_col = None
        s.sns_id_col = None
        s.sns_gene_col = None
        s.hxt_score_col = None
        s.hxt_id_col = None
        s.hxt_gene_col = None

        #This will be a dictionary of species. Each species is a list of tuples
        # of (gene, sensor_score, transporter_score). Dictionary entry will
        # be the specified ID field with spaces replaced with underscores.
        # Also converted to all lowercase.
        #TODO In retrospect, there is no reason to use tuples here. I have
        # to modify the tuple anyway in the transporter portion, so I should
        # use a list since they are much faster anyway. Consider it.
        s.scores = {}
        s.scores_work = {}

        #If I've got to open any more files, make this a function
        #Also note that it's possible that turning the csv reader into a 
        # list of lists is significantly slower than just iterating over
        # the CSV reader object. Keep in mind if performance issues arise.
        #TODO Sniffer for delimiters, quotechar, etc.
        try:
            s.sns_file = open(sns_filename, 'r')
        except IOError:
            print('Cannot find/open sensor HMMER file')
        else:
            s.sns_csvrdr = csv.reader(s.sns_file, delimiter=s.delimiter, 
                                      quotechar=s.quotechar)
            s.sns_csvlist = list(s.sns_csvrdr)
            #Removing label entry. TODO Implement sniffer usage later, it can 
            # do this kind of functionality in a much nicer way.
            del s.sns_csvlist[0]

        try:
            s.hxt_file = open(hxt_filename, 'r')
        except IOError:
            print('Cannot find/open transporter HMMER file')
        else:
            s.hxt_csvrdr = csv.reader(s.hxt_file, delimiter=s.delimiter, 
                                      quotechar=s.quotechar)
            s.hxt_csvlist = list(s.hxt_csvrdr)
            #TODO Sniffer here too
            del s.hxt_csvlist[0]

    def setCSVVars(s, sns_score_col, sns_id_col, sns_gene_col,
                   hxt_score_col, hxt_id_col, hxt_gene_col):
        #Input sanitization will not happen here. Instead, we will
        # catch the TypeError exception in the for-loop. Ask forgivenss
        # not permission.

        #No way to make sure col values is in bounds of CSV without going 
        # through entire file once before which is slow and pointless. Throw 
        # exceptions for OOB errors at actual operation time. Forgivenss not
        # permission.  !! Actually this is now doable due to conversion of
        # CSV reader into a list. Decide best method for this.

        s.sns_score_col = sns_score_col
        s.sns_id_col = sns_id_col
        s.sns_gene_col = sns_gene_col
        s.hxt_score_col = hxt_score_col
        s.hxt_id_col = hxt_id_col
        s.hxt_gene_col = hxt_gene_col

    def pairScores(s, truncate_ident=0):
        try:
            for entry in s.sns_csvlist:
                ident = entry[s.sns_id_col].replace(' ', '_').lower()
                if truncate_ident:
                    ident_split = ident.split('_')
                    ident = ''
                    for x in range(truncate_ident):
                        ident += ident_split[x] + '_'
                gene = entry[s.sns_gene_col].replace(' ', '_').lower()
                score = entry[s.sns_score_col]
                if ident not in s.scores:
                    s.scores[ident] = [(gene, score,)]
                #Check to make sure it's not a dupe
                #Do we even care if there are dupes?
                elif(gene, score,) not in s.scores[ident]:
                    s.scores[ident].append((gene, score,))
                #It's a dupe
                else:
                    pass
        except ValueError:
            print('CSV column values for sensor file must be integers.')
        except TypeError:
            print('Remember to load column values using setCSVVars(...)!')
        except IndexError:
            print('CSV column values for sensor file must be between 0 and the'
                  ' total number of columns-1, inclusive.  The provided file' 
                  ' has ', len(s.sns_csvlist[0]), ' columns total.')

        try:
            for entry in s.hxt_csvlist:
                ident = entry[s.hxt_id_col].replace(' ', '_').lower()
                if truncate_ident:
                    ident_split = ident.split('_')
                    ident = ''
                    for x in range(truncate_ident):
                        ident += ident_split[x] + '_'
                gene = entry[s.hxt_gene_col].replace(' ', '_').lower()
                score = entry[s.hxt_score_col]
                #Will be an entire unpaired species in this case
                if ident not in s.scores:
                    s.scores[ident] = [(gene, None, score,)]
                #Check to make sure it's not a dupe
                else:
                    #This list comp must return something then
                    gs_index = [x for x, genescore in 
                                enumerate(s.scores[ident])
                                if genescore[0] == gene]
                    #Gene doesn't exist for this species in sensor list
                    if len(gs_index) == 0:
                        s.scores[ident].append((gene, None, score,))
                    #It's a dupe (already been appended)
                    elif len(s.scores[ident][gs_index[0]]) > 2:
                        pass
                    #It's not, go ahead and append
                    else:
                        # When it becomes a list, do this instead
                        #s.scores[ident][gs_index].append(score)
                        s.scores[ident][gs_index[0]] += (score,)
            #Need to append 'None' to all the genes found in sensor list
            # but not in hxt list. No way in hell is the cleanest solution.
            # This in list comp form would be gross looking.
            for i in s.scores:
                for g in enumerate(s.scores[i]):
                    if len(s.scores[i][g[0]]) == 2:
                        s.scores[i][g[0]] += (None,)
        except ValueError:
            print('CSV column values for transporter file must be integers.')
        except IndexError:
            print('CSV column values for transporter file must be between 0'
                  ' and the total number of columns-1, inclusive.  The' 
                  ' provided file has ', len(s.hxt_csvlist[0]), ' columns'
                  ' total.')

    #Will need to make this take some sort of temporary 'acting' scores dict
    # that can be modified with member functions.
    def makePlots(s):
        pass

    #This member function will detect whether to do single species or multiple
    # species based on the name. Perhaps have a member generate full genus
    # entries as seperate entities in the dict?
    def showPlot(s, name):
        x = [i[1] for i in s.scores[name]]
        y = [i[2] for i in s.scores[name]]

        x = [0 if i==None else i for i in x]
        y = [0 if i==None else i for i in y] 

        plt.scatter(x, y)
        plt.plot([0, 900], [0, 900], 'k-')
        plt.show()

    #This operates on the working scores array. 
    # string=(<string>, ...): A tuple of the strings to search for
    # action='remove', 'keep': Whether to remove all entries in which
    #       string appears in field or remove ALL BUT the entries in which
    #       string appears in field
    # array='working', 'base': Whether to use the base scores array or the 
    #       working scores array as the resource to operate on. In, either 
    #       case, the old working scores array is overwritten.
    #
    #By default, this function removes all entries from the base set to add
    # to the working scores array, effectively clearing the working scores 
    # array.
    #TODO: Implement a field parameter to be able to search fields other
    #       than just species identifier. Gene name?
    def filterPairsByString(s, strings=('',),  action='remove', array='base'):
        temp_scores = {}
        
        if array=='base': source = s.scores
        elif array=='working': source = s.scores_work
        else: 
            print('The array parameter must either be \'base\' or'
                    ' \'working\'.')
            return

        for ident in source:
            for string in strings:
                if string in ident:
                    if action=='keep':
                        temp_scores[ident] = copy.deepcopy(source[ident]) 
                    elif action!='remove':
                        print('The action parameter must be \'keep\' or'
                              ' \'remove\'.')
                        return
                elif action=='remove':
                    temp_scores[ident] = copy.deepcopy(source[ident])

        s.scores_work = temp_scores

