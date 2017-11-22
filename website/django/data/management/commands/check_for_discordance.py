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
        consistency_list_positive = ['pathogenic', 'likely pathogenic', '+', '+?']
        consistency_list_negative = ['benign', 'likely benign', 'unclassified', 'uncertain', '-', '-?', '?', '.']
        return self._determine_consistency_or_concordance(consistency_list_positive, consistency_list_negative, output_row, 'consistency')


    def _determine_discordance(self, output_row):
        discordance_list_positive = ['pathogenic', 'likely pathogenic', '+', '+?']
        discordance_list_negative = ['benign', 'likely benign', '-', '-?']
        clinvar_concordance, lovd_concordance, clinvar_and_lovd_concordance = self._determine_consistency_or_concordance(discordance_list_positive, discordance_list_negative, output_row, 'concordance')
        # return discordance, not concordance (opposite)
        clinvar_discordance = False if clinvar_concordance is None else not clinvar_concordance
        lovd_discordance = False if lovd_concordance is None else not lovd_concordance
        clinvar_and_lovd_discordance = False if clinvar_and_lovd_concordance is None else not clinvar_and_lovd_concordance
        return clinvar_discordance, lovd_discordance, clinvar_and_lovd_discordance


    def _determine_consistency_or_concordance(self, list_positive, list_negative, output_row, type):
        clinvar_positive = False
        lovd_positive = False
        clinvar_negative = False
        lovd_negative = False
        clinvar_consistency_or_concordance = None
        lovd_consistency_or_concordance = None
        clinvar_and_lovd_consistency_or_concordance = None

        if not is_empty(output_row["Clinical_Significance_ClinVar"]):
            for clinvar_sig in output_row["Clinical_Significance_ClinVar"].split(','):
                if clinvar_sig in list_positive:
                    clinvar_positive = True
                elif clinvar_sig in list_negative:
                    clinvar_negative = True
        
        if not is_empty(output_row["Variant_effect_LOVD"]):
            for lovd_sig in output_row["Variant_effect_LOVD"].split(','):
                lovd_sig = lovd_sig.split('/')[1]
                if lovd_sig in list_positive:
                    lovd_positive = True
                elif lovd_sig in list_negative:
                    lovd_negative = True


        # Determine clinvar consistency/concordance

        if clinvar_positive is True and clinvar_negative is True:
            clinvar_consistency_or_concordance = False
        elif type == "concordance":
            if clinvar_positive is True and clinvar_negative is False:
                clinvar_consistency_or_concordance = True
            elif clinvar_negative is True and clinvar_positive is False:
                clinvar_consistency_or_concordance = True
        elif type == "consistency":
            # clinvar consistency only cares if all reports are positive (consistently actionable)
            if clinvar_positive is True and clinvar_negative is False:
                clinvar_consistency_or_concordance = True
            else:
                clinvar_consistency_or_concordance = False
        else:
            print "Must check either concordance or consistency."
            raise Exception
            sys.exit(1)


        # Determine lovd consistency/concordance

        if lovd_positive is True and lovd_negative is True:
            lovd_consistency_or_concordance = False
        elif lovd_positive is True and lovd_negative is False:
            lovd_consistency_or_concordance = True
        elif lovd_negative is True and lovd_positive is False:
            lovd_consistency_or_concordance = True

        
        # Determine combined consistency/concordance

        if clinvar_positive is True and lovd_negative is True:
            clinvar_and_lovd_consistency_or_concordance = False
        elif clinvar_negative is True and lovd_positive is True:
            clinvar_and_lovd_consistency_or_concordance = False
        elif clinvar_positive is True and clinvar_negative is True:
            clinvar_and_lovd_consistency_or_concordance = False
        elif lovd_positive is True and lovd_negative is True:
            clinvar_and_lovd_consistency_or_concordance = False
        elif clinvar_positive is True and lovd_positive is True:
            clinvar_and_lovd_consistency_or_concordance = True
        elif clinvar_negative is True and lovd_negative is True:
            clinvar_and_lovd_consistency_or_concordance = True
        else:
            vals = [clinvar_positive, clinvar_negative, lovd_positive, lovd_negative]
            num_true = 0
            for val in vals:
                if val is True:
                    num_true += 1
            if num_true <= 1:
                clinvar_and_lovd_consistency_or_concordance = True

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

            clinvar_actionable, lovd_consistency, clinvar_and_lovd_consistency = self._determine_consistency(output_row)
            output_row['Actionable_ClinVar'] = clinvar_actionable
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
            'file': '/tmp/discordance.tsv',
            'class': CurrentVariant,
            'variant_fields': ['Genomic_Coordinate_hg38', 'HGVS_cDNA',
                       'Allele_frequency_ExAC', 'Allele_frequency_1000_Genomes'],
            'fields_of_interest': ['Genomic_Coordinate_hg38', 'HGVS_cDNA', 'Submitter_ClinVar',
                       'Clinical_Significance_ClinVar', 'Actionable_ClinVar',
                       'Discordance_ClinVar', 'Submitters_LOVD', 'Variant_effect_LOVD',
                       'Consistency_LOVD', 'Discordance_LOVD', 'Consistency_LOVD_And_ClinVar',
                       'Discordance_LOVD_And_ClinVar', 'Allele_frequency_ExAC',
                       'Allele_frequency_1000_Genomes']
        }

        f = open(meta['file'], 'w+')
        writer = csv.writer(f, encoding='utf-8', delimiter='\t')
        writer.writerow( meta['fields_of_interest'] )
        for obj in meta['class'].objects.all():
            obj_data = self._collect_data_for_variant(obj, meta)
            if not obj_data:
                continue
            else:
                print obj_data
                print '\n'
                row = []
                for field in meta['fields_of_interest']:
                    row.append(obj_data[field])
                writer.writerow(row)
        f.close()
        print 'Data written to %s' % meta['file']
