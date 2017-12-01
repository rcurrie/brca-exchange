# coding=utf-8

"""
Determines consistency and discordance among data in database for submissions from ClinVar and LOVD.
"""

from django.core.management.base import BaseCommand, CommandError
from django.db import connection
from data.models import CurrentVariant, Report, Variant
from django.db import transaction
import unicodecsv as csv


EMPTY = '-'


def is_empty(value):
    return value == EMPTY


class Command(BaseCommand):
    help = 'Determine discordance in data and output to a tsv file.'


    def _append_value(self, obj, field, value):
        if value is None:
            value = EMPTY
        if is_empty(obj[field]):
            obj[field] = value
        else:
            obj[field] += ', ' + value
        return obj


    def _determine_consistency(self, output_row):
        consistency_list_positive = ['pathogenic', 'likely_pathogenic', 'probable-pathogenic', '+', '+?']
        consistency_list_negative = ['benign', 'likely_benign', 'unclassified', 'uncertain_significance', 'probably_not_pathogenic',
                                     'no_known_pathogenicity', 'not_provided', 'variant_of_unknown_significance', '-', '-?', '?', '.']
        return self._determine_consistency_or_concordance(consistency_list_positive, consistency_list_negative, output_row)


    def _determine_discordance(self, output_row):
        discordance_list_positive = ['pathogenic', 'likely_pathogenic', 'probable-pathogenic', '+', '+?']
        discordance_list_negative = ['benign', 'likely_benign', 'probably_not_pathogenic', '-', '-?']
        clinvar_concordance, lovd_concordance, clinvar_and_lovd_concordance = self._determine_consistency_or_concordance(discordance_list_positive, discordance_list_negative, output_row)
        # return discordance, not concordance (opposite)
        clinvar_discordance = False if clinvar_concordance is None else not clinvar_concordance
        lovd_discordance = False if lovd_concordance is None else not lovd_concordance
        clinvar_and_lovd_discordance = False if clinvar_and_lovd_concordance is None else not clinvar_and_lovd_concordance
        return clinvar_discordance, lovd_discordance, clinvar_and_lovd_discordance


    def _determine_clinvar_or_lovd_consistency_or_concordance(self, positive, negative):
        if positive is True and negative is True:
            return False
        elif positive is True and negative is False:
            return True
        elif negative is True and positive is False:
            return True
        else:
            return None


    def _determine_clinvar_and_lovd_consistency_or_concordance(self, clinvar_positive, clinvar_negative, lovd_positive, lovd_negative):
        if clinvar_positive is True and lovd_positive is True and clinvar_negative is False and lovd_negative is False:
            return True
        elif clinvar_negative is True and lovd_negative is True and clinvar_positive is False and lovd_positive is False:
            return True
        else:
            if [clinvar_positive, clinvar_negative, lovd_positive, lovd_negative].count(True) <= 1:
                return True
            else:
                return False


    def _determine_consistency_or_concordance(self, list_positive, list_negative, output_row):
        clinvar_positive = False
        lovd_positive = False
        clinvar_negative = False
        lovd_negative = False

        if not is_empty(output_row["Clinical_Significance_ClinVar"]):
            for clinvar_sig in output_row["Clinical_Significance_ClinVar"].split(','):
                clinvar_sig = clinvar_sig.strip().lower()
                if clinvar_sig in list_positive:
                    clinvar_positive = True
                elif clinvar_sig in list_negative:
                    clinvar_negative = True
        
        if not is_empty(output_row["Variant_effect_LOVD"]):
            for lovd_sig in output_row["Variant_effect_LOVD"].split(','):
                lovd_sig = lovd_sig.strip().lower()
                lovd_sig = lovd_sig.split('/')[1]
                if lovd_sig in list_positive:
                    lovd_positive = True
                elif lovd_sig in list_negative:
                    lovd_negative = True

        clinvar_consistency_or_concordance = self._determine_clinvar_or_lovd_consistency_or_concordance(clinvar_positive, clinvar_negative)
        lovd_consistency_or_concordance = self._determine_clinvar_or_lovd_consistency_or_concordance(lovd_positive, lovd_negative)
        clinvar_and_lovd_consistency_or_concordance = self._determine_clinvar_and_lovd_consistency_or_concordance(clinvar_positive, clinvar_negative, lovd_positive, lovd_negative)

        return clinvar_consistency_or_concordance, lovd_consistency_or_concordance, clinvar_and_lovd_consistency_or_concordance


    def _collect_data_for_variant(self, obj, meta):
        # get reports related to variant
        variant = Variant.objects.get(id=obj.id)
        reports = variant.report_set.all()
        
        # build dictionary of fields of interest
        output_row = {}
        for field_to_write in meta['fields_of_interest']:
            output_row[field_to_write] = EMPTY    
        
        # determine if variant is of interest
        # TODO: ignore cases where clinvar and lovd have data from same submitter
        relevant_reports = []
        for report in reports:
            if report.Source == "ClinVar":
                relevant_reports.append(report)
            elif report.Source == "LOVD":
                relevant_reports.append(report)

        # if variant is of interest, fill out fields of interest
        if len(relevant_reports) >= 2:

            for field in meta['variant_fields']:
                output_row[field] = getattr(variant, field)
            for report in relevant_reports:
                if report.Source == "ClinVar":
                    significance = report.Clinical_Significance_ClinVar.lower()
                    output_row = self._append_value(output_row, 'Submitter_ClinVar', report.Submitter_ClinVar)
                    output_row = self._append_value(output_row, 'Clinical_Significance_ClinVar', significance)
                elif report.Source == "LOVD":
                    significance = report.Variant_effect_LOVD
                    output_row = self._append_value(output_row, 'Submitters_LOVD', report.Submitters_LOVD)
                    output_row = self._append_value(output_row, 'Variant_effect_LOVD', significance)

            clinvar_consistency, lovd_consistency, clinvar_and_lovd_consistency = self._determine_consistency(output_row)
            output_row['Consistency_ClinVar'] = clinvar_consistency
            output_row['Consistency_LOVD'] = lovd_consistency
            output_row["Consistency_LOVD_And_ClinVar"] = clinvar_and_lovd_consistency

            clinvar_discordance, lovd_discordance, clinvar_and_lovd_discordance = self._determine_discordance(output_row)
            output_row['Discordance_ClinVar'] = clinvar_discordance
            output_row['Discordance_LOVD'] = lovd_discordance
            output_row["Discordance_LOVD_And_ClinVar"] = clinvar_and_lovd_discordance

            return output_row
        else:
            return False


    def handle(self, *args, **options):
        meta = {
            'file': '/Users/zfisch/UCSC/data/discordance/discordance.tsv',
            'variant_fields': ['Genomic_Coordinate_hg38', 'HGVS_cDNA',
                       'Allele_frequency_ExAC', 'Allele_frequency_1000_Genomes'],
            'fields_of_interest': ['Genomic_Coordinate_hg38', 'HGVS_cDNA', 'Submitter_ClinVar',
                       'Clinical_Significance_ClinVar', 'Consistency_ClinVar',
                       'Discordance_ClinVar', 'Submitters_LOVD', 'Variant_effect_LOVD',
                       'Consistency_LOVD', 'Discordance_LOVD', 'Consistency_LOVD_And_ClinVar',
                       'Discordance_LOVD_And_ClinVar', 'Allele_frequency_ExAC',
                       'Allele_frequency_1000_Genomes']
        }

        f = open(meta['file'], 'w+')
        writer = csv.writer(f, encoding='utf-8', delimiter='\t')
        writer.writerow( meta['fields_of_interest'] )
        print 'Building discordance.tsv -- this takes a while.'
        for obj in CurrentVariant.objects.all().exclude(Change_Type__name='deleted'):
            obj_data = self._collect_data_for_variant(obj, meta)
            if not obj_data:
                continue
            else:
                row = []
                for field in meta['fields_of_interest']:
                    row.append(obj_data[field])
                writer.writerow(row)
        f.close()
        print 'Data written to %s' % meta['file']
