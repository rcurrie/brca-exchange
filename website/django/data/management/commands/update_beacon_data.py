from django.core.management.base import BaseCommand, CommandError
from django.db import connection, transaction
from data.models import Variant, CurrentVariant, DataRelease, ChangeType, Report
import logging
import time
import json
import requests
from django.core.management import call_command
import pdb

URL_MAX = 1000  # The maximum allowed length of HTTP request URLs
BEACON_RESPONSE_URL = "https://beacon-network.org/api/responses" # URL for variant requests in the beacon network
BEACON_BEACONS_URL = "https://beacon-network.org/api/beacons" # URL for beacon list requests
BRCA_GA4GH_URL = "http://brcaexchange.cloudapp.net/backend/data/ga4gh/v0.6.0a7/" # URL for BRCA Exchange GA4GH instance

class Command(BaseCommand):
    help = 'Submits queries against the Beacon Network for all the single nucleotide variants and adds results to the database.'

    def add_arguments(self, parser):
        parser.add_argument('-t', '--tries', default=5, type=int,
                            help="The number of times to attempt each beacon network request")
        parser.add_argument('-p', '--test', default=False,
                            help="Work in a 10 base range")
        # parser.add_argument('-o', '--output', type=FileType('w'), help='Output file to store results from queries.')

    def handle(self, *args, **options):
        # output_file = options['output']
        tries = options['tries']
        test = options['test']
        
        # Get latest variants
        variants = CurrentVariant.objects.all().exclude(Change_Type__name='deleted')

        # TODO: make function
        # Get current list of beacons
        for i in range(tries):
            try:
                beacons = json.loads(requests.get(BEACON_BEACONS_URL).content)
                beacon_ids = map(lambda x: str(x["id"]), beacons)
                for beacon in beacons:
                    if beacon['aggregator'] == True:
                        beacon_ids.remove(str(beacon['id']))
                beacon_id_str = "[" + ",".join(beacon_ids) + "]"
                break
            except:
                logging.error("Failed to access the Beacon Network")
                if i+1 == tries:
                    raise
            time.sleep(20)
            
        # TODO: make functions
        for variant in variants:
            if test:
                BRCA1_CHR = "chr17" 
                BRCA1_START = 41223107
                BRCA1_END = 41223109

                BRCA2_CHR = "chr13"
                BRCA2_START = 32929255
                BRCA2_END = 32929257

                GENE_BUFFER = 5

                if variant.Chr == "17":
                    if (variant.Hg37_Start < (BRCA1_START - GENE_BUFFER)) or (variant.Hg37_End > (BRCA1_END + GENE_BUFFER)):
                        continue
                if variant.Chr == "13":
                    if (variant.Hg37_Start < (BRCA2_START - GENE_BUFFER)) or (variant.Hg37_End > (BRCA2_END + GENE_BUFFER)):
                        continue
            # beacon_data = query_beacons_for_variant(variant)

            # @transaction.atomic

            if len(variant.Ref) != 1 or len(variant.Alt) != 1:
                save_beacon_data_for_skipped_variant(variant, "Variant is not a SNP")
                logging.info("Variant {} is not a SNP: skipping".format(variant.id))
                continue
            

            # TODO: ensure -1 is needed for start
            request = "?chrom={}&pos={}&allele={}&ref=GRCh37".format(
                          variant.Chr, variant.Hg37_Start-1,
                          variant.Alt)


            if len(request) > URL_MAX:
                save_beacon_data_for_skipped_variant(variant, "Variant URL is too long")
                logging.info("Skipping variant {}: URL too long".format(variant.id))
                continue

            variant_dictionary = query_beacons_one_request(beacons, variant.id, request, tries)
            website_url = "https://beacon-network.org/#/search?chrom={}&pos={}&ref={}&allele={}&rs=GRCh37".format(
                              variant.Chr, variant.Hg37_Start,
                              variant.Ref, variant.Alt)
            variant_dictionary["url"] = website_url

            print variant_dictionary


#     # update materialized view of current variants
#     with connection.cursor() as cursor:
#         cursor.execute("REFRESH MATERIALIZED VIEW currentvariant")


def query_beacons_one_request(beacons, variant, request, tries):
    """
    Queries the given list of beacons for the given variant in a single request
    """
    contents = ''
    for i in range(tries):
        try:
            r = requests.get(BEACON_RESPONSE_URL + request, timeout=180)
            logging.debug("Requesting {}{}".format(BEACON_RESPONSE_URL, request))
            contents = json.loads(r.content)
            break
        except:
            # Report failure if we've run out of tries
            if i+1 == tries:
                logging.error("Failed request: variant: {}".format(variant))

    return build_variant_dictionary(variant, contents)


# class QueryBeacons (threading.Thread):
#     """
#     Starts a thread that queries each of the beacons for the given variant, skipping aggregators
#     """

#     def __init__(self, beacons, variant, request, tries):
#         threading.Thread.__init__(self)
#         self.beacons = beacons
#         self.variant = variant
#         self.request = request
#         self.tries = tries

#     def run(self):

#         query_beacons(self.beacons, self.variant, self.request, self.tries)


# def query_beacons(beacons, variant, request, tries):
#     """
#     Queries the given list of beacons for the given variant one beacon at a time, skipping
#     the beacon-of-beacons. 
#     """
#     beacon_responses = []
#     for beacon in beacons:
#     logging.info("Requesting variant from beacon: {}".format(beacon["id"]))        
#         if beacon["id"] == "bob":
#             continue

#         for i in range(tries):
#             try:
#                 r = requests.get(BEACON_RESPONSE_URL + "/" +  beacon["id"] + request, timeout=20)
#                 logging.debug("Requesting {}/{}{}".format(BEACON_RESPONSE_URL, beacon["id"], request))
#                 content = json.loads(r.content)
#                 beacon_responses.append(content)
#                 break

#             except:
#                 # Report failure if we've run out of tries
#                 if i+1 == tries:
#                     logging.error("Failed request: variant: {}, beacon: {}".format(variant, beacon["id"]))

    build_variant_dictionary(variant, beacon_responses)


def build_variant_dictionary(variant, beacon_responses):
    """
    Adds data from beacon_responses to the variant dictionary entry for variant
    """
    variant_dictionary = {"beaconNames": [], "beaconIDs": [], "numBeacons": 0, "notFound": 0, "notApplicable": 0, "error": None}
    for response in beacon_responses:
        if 'response' in response:
            if response['response']:
                if response['beacon']['aggregator']:
                    continue

                variant_dictionary['beaconNames'].append(response['beacon']['name'])
                variant_dictionary['beaconIDs'].append(response['beacon']['id'])
                variant_dictionary['numBeacons'] += 1
            elif response['response'] == None:
                variant_dictionary['notApplicable'] += 1
            elif response['response'] == False:
                variant_dictionary['notFound'] += 1
        else:
            variant_dictionary['notApplicable'] += 1
            logging.debug("No response from {}: {}".format(response['beacon']["id"], str(response)))
    return variant_dictionary


def save_beacon_data_for_skipped_variant(variant, error):
    variant_dictionary = {"beaconNames": [], "beaconIDs": [], "numBeacons": 0, "notFound": 0, "notApplicable": 0, "url": None, "error": error}
    print variant_dictionary