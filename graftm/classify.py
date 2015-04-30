import json

class Classify:
    def __init__(self,taxonomy):
        self.taxonomy=self.readRefpkgTax(taxonomy)

    def readRefpkgTax(self, taxonomy_file):
        ## Read in the taxonomic description of the tree within the refpkg
        taxonomy_hash={}
        for line in [[y for y in x.rstrip().split(',') if y] for x in open(taxonomy_file).readlines() if not x.startswith('tax_id')]:
            taxonomy_hash[line[0]] = ['Root']+line[5:]
        return taxonomy_hash

    def assignPlacement(self, placement_json_path, cutoff, type):
        ## Function that reads in classification and returns a 'guppy classify'
        ## like file
        all_placements_reads={}

        def reduceTaxString(placement_hashes, threshold):
            confidences=[x['c'] for x in placement_hashes]
            tax_that_meets_threshold={'placement':[],
                                      'confidence':[]}
            for i in range(0,len(max([x['p'] for x in placement_hashes]))):
                cumil_confidence={}
                try:
                    for idx, item in enumerate([x['p'][i] for x in placement_hashes]):
                        if item in cumil_confidence.keys():
                            cumil_confidence[item]+=confidences[idx]
                        else:
                            cumil_confidence[item]=confidences[idx]
                    best_place=max([(value, key) for key, value in cumil_confidence.items()])
                except IndexError:
                    continue
                if best_place[0]>threshold:
                    tax_that_meets_threshold['placement'].append(best_place[1])
                    tax_that_meets_threshold['confidence'].append(best_place[0])
                else:
                    return tax_that_meets_threshold
            if tax_that_meets_threshold['placement']:
                return tax_that_meets_threshold
            else:
                raise Exception("Programming error.")

        def consolidatePlacements(placement_list, cutoff, lwr_idx, c_idx):
            seen={}
            for placement in placement_list:
                rank=placement[0]
                confidence=placement[lwr_idx]
                if placement[0] not in seen:
                    taxonomy_string=self.taxonomy[rank]
                    seen[rank]={'c':confidence,
                                'p':taxonomy_string}
                else:
                    seen[rank]['c']+=confidence
            if len(seen)==1 and seen.items()[0][1]['c']>=0.98: # If there is one entry in seen, and that entry has full confidence.
                return {'placement': seen.items()[0][1]['p'], 'confidence': [seen.items()[0][1]['c']]*len(seen.items()[0][1]['p'])} # Return that tax            
            elif len(seen)>1: # If there is more than one entry
                return reduceTaxString(seen.values(), cutoff)
            else:
                raise Exception("Programming Error: Classify; assignPlacement; consolidatePlacements")

        placement_hash=json.load(open(placement_json_path)) # read in placement json
        try: # Search for the idx of the like field ratio and classification 
            lwr_idx=placement_hash['fields'].index('like_weight_ratio')
            c_idx=placement_hash['fields'].index('classification')
        except ValueError: # If they can't be found, FAAAAAIIILL.
            raise Exception('Fatal error in refpkg, classification or like_weight_ratio fields missing')

        for placement_group in placement_hash['placements']: # for each placement
            best_place=consolidatePlacements(placement_group['p'], cutoff, lwr_idx, c_idx) # Find the best placement      
            if best_place: # if it exists
                reads=[x[0] for x in placement_group['nm']] # make a list of the reads assigned to that placement
                for read in reads: # and for each read
                    if read[-1] in all_placements_reads.keys(): # Sort each read by its file index and enter it into hash, with the best placement as the value
                        all_placements_reads[read[-1]][read[:-2]] = best_place 
                    else:
                        all_placements_reads[read[-1]] = {read[:-2]: best_place}
            else: # If the best placement doesn't exist, fail.
                raise Exception("Programming Error: Failed to retrieve taxonomy for %s"  % (' '.join([x[0] for x in placement_group['nm']])))

        return all_placements_reads
