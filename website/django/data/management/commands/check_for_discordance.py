# coding=utf-8

from django.core.management.base import BaseCommand, CommandError
from django.db import connection
from data.models import CurrentVariant, Report, Variant
from django.db import transaction
import unicodecsv as csv
import pdb

# - the variant (HGVS cDNA coordinates would probably be easiest for most people)
# - a list of the submitters who've submitted this variant to ClinVar
# - clinical significance according to each of these submitters (if this field and the last field are comma-delimited with parallel indices, such as the Nth submitter corresponds to the Nth clinical significance entry, that should work as a general solution)
# - a binary flag indicating whether or not the variant is consistently assessed as actionable (pathogenic or likely pathogenic) according to this list of assessments
# - a binary flag indicating if the list of assessments is truly discordant (i.e. the variant has been classified as both (pathogenic, likely pathogenic) and (benign, likely benign)).
# - a list of the submitters who've submitted this variant to LOVD
# - the curated variant effect, according to each of these submissions (e.g. given a variant effect such as '?/.', the curator's assessment is the second one, '.').
# - whether or not those variant effect reports are consistent.  More on this below.
# - whether or not those variant effect reports are truly discordant.  More on this below
# - Whether the variant is assessed consistently across everything in ClinVar AND LOVD
# - Whether the variant has any truly discordant assessments in ClinVar and/or LOVD, with all assessments taken together
# - overall allele frequency from ExAC, if available
# - overall allele frequency from 1000 Genomes, if available


EMPTY = '-'


def is_empty(value):
    return value == EMPTY


class Command(BaseCommand):
    help = 'Determine discordance in data and output to a tsv file.'

    """
    NOTE: Let's look at this for variants that have two or more distinct
    submissions to ClinVar and LOVD collectively. One nuance is that if
    the same submitter submits the same variant to ClinVar and LOVD, it's
    not two distinct submissions but one submission reported twice. 
    I know that ENIGMA has some submissions in both LOVD and ClinVar,
    and there are probably other groups that do as well."
    """

    def _append_value(self, obj, field, value):
        if is_empty(obj[field]):
            obj[field] = value
        else:
            obj[field] += ', ' + value
        return obj


    def _update_clinvar_discordance(self, discordance, significance, first_clinvar_significance):
        if discordance is True:
            return True, first_clinvar_significance
        elif first_clinvar_significance is None:
            discordance = False
            first_clinvar_significance = significance
        elif "pathogenic" in first_clinvar_significance:
            if "benign" in significance:
                discordance = True
        elif "benign" in first_clinvar_significance:
            if "pathogenic" in significance:
                discordance = True
        else:
            # if neither pathogenic or benign (e.g. unknown), update first_significance
            first_clinvar_significance = significance
        return discordance, first_clinvar_significance


    def _update_clinvar_actionable(self, significance, actionable):
        if actionable is False:
            return False
        if "pathogenic" in significance:
            actionable = True
        else:
            actionable = False
        return actionable


    def _update_lovd_consistency(self, consistency, significance, first_lovd_significance):
        # Curator's assessment is the value after the / (e.g. '?/.' would be '.')
        consistency_list_positive = ['pathogenic', 'likely pathogenic', '+', '+?']
        consistency_list_negative = ['benign', 'likely benign', 'unclassified', 'uncertain', '-', '-?', '?', '.']
        
        significance = significance.split('/')[1].lower()
        
        if consistency is False:
            return False, first_lovd_significance
        elif consistency is None:
            if significance in consistency_list_negative or significance in consistency_list_positive:
                consistency = True
                first_lovd_significance = significance
            else:
                consistency = None
                first_lovd_significance = None
        elif first_lovd_significance in consistency_list_positive:
            if significance in consistency_list_negative:
                consistency = False
        else:
            if significance in consistency_list_positive:
                consistency = False
        return consistency, first_lovd_significance

    def _update_lovd_discordance(self, discordance, significance, first_lovd_significance):
        discordance_list_positive = ['pathogenic', 'likely pathogenic', '+', '+?']
        discordance_list_negative = ['benign', 'likely benign', '-', '-?']

        significance = significance.split('/')[1].lower()
        
        if discordance is True:
            return True, first_lovd_significance
        elif first_lovd_significance is None:
            if significance in discordance_list_positive or significance in discordance_list_negative:
                discordance = False
                first_lovd_significance = significance
            else:
                discordance = None
                first_lovd_significance = None
        elif first_lovd_significance in discordance_list_positive:
            if significance in discordance_list_negative:
                discordance = True
        elif first_lovd_significance in discordance_list_negative:
            if significance in discordance_list_positive:
                discordance = True
        return discordance, first_lovd_significance


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
        relevant_reports_count = 0
        for report in reports:
            if report.Source == "ClinVar":
                relevant_reports_count += 1
            elif report.Source == "LOVD":
                relevant_reports_count += 1

        # if variant is of interest, fill out fields of interest
        if relevant_reports_count >= 2:
            clinvar_discordance = None
            clinvar_actionable = None
            first_clinvar_significance = None
            lovd_consistency = None
            lovd_discordance = None
            first_lovd_significance = None
            for field in meta['variant_fields']:
                output_row[field] = getattr(variant, field)
            for report in reports:
                if report.Source == "ClinVar":
                    significance = report.Clinical_Significance_ClinVar.lower()
                    clinvar_actionable = self._update_clinvar_actionable(significance, clinvar_actionable)
                    first_clinvar_significance = None
                    clinvar_discordance, first_clinvar_significance = self._update_clinvar_discordance(clinvar_discordance, significance, first_clinvar_significance)
                    output_row = self._append_value(output_row, 'Submitter_ClinVar', report.Submitter_ClinVar)
                    output_row = self._append_value(output_row, 'Clinical_Significance_ClinVar', significance)
                elif report.Source == "LOVD":
                    significance = report.Variant_effect_LOVD
                    lovd_consistency, first_lovd_significance = self._update_lovd_consistency(lovd_consistency, significance, first_lovd_significance)
                    first_lovd_significance = None
                    lovd_discordance, first_lovd_significance = self._update_lovd_discordance(lovd_discordance, significance, first_lovd_significance)
                    output_row = self._append_value(output_row, 'Submitters_LOVD', report.Submitters_LOVD)
                    output_row = self._append_value(output_row, 'Variant_effect_LOVD', significance)
            output_row['Consistency_LOVD'] = lovd_consistency
            output_row['Discordance_LOVD'] = lovd_discordance
            output_row['Actionable_ClinVar'] = clinvar_actionable
            output_row['Discordance_ClinVar'] = clinvar_discordance

            consistency_list_positive = ['pathogenic', 'likely pathogenic', '+', '+?']
            consistency_list_negative = ['benign', 'likely benign', 'unclassified', 'uncertain', '-', '-?', '?', '.']

            if lovd_consistency is False:
                output_row["Consistency_LOVD_And_ClinVar"] = False
            else:
                clinvar_positive = None
                clinvar_negative = None
                lovd_positive = None
                lovd_negative = None
                for clinvar_sig in output_row["Clinical_Significance_ClinVar"].split(','):
                    if clinvar_sig in consistency_list_positive:
                        clinvar_positive = True
                    elif clinvar_sig in consistency_list_negative:
                        clinvar_negative = True
                for lovd_sig in output_row["Variant_effect_LOVD"].split(','):
                    if lovd_sig in consistency_list_positive:
                        lovd_positive = True
                    elif lovd_sig in consistency_list_negative:
                        lovd_negative = True
                if clinvar_positive is True and clinvar_negative is True:
                    output_row["Consistency_LOVD_And_ClinVar"] = False
                elif clinvar_positive is True and lovd_negative is True:
                    output_row["Consistency_LOVD_And_ClinVar"] = False
                elif lovd_positive is True and clinvar_negative is True:
                    output_row["Consistency_LOVD_And_ClinVar"] = False
                elif lovd_positive is True and lovd_negative is True:
                    output_row["Consistency_LOVD_And_ClinVar"] = False
                else:
                    output_row["Consistency_LOVD_And_ClinVar"] = True

            discordance_list_positive = ['pathogenic', 'likely pathogenic', '+', '+?']
            discordance_list_negative = ['benign', 'likely benign', '-', '-?']

            if lovd_discordance is True:
                output_row["Discordance_LOVD_And_ClinVar"] = True
            elif clinvar_discordance is True:
                output_row["Discordance_LOVD_And_ClinVar"] = True
            else:
                clinvar_positive = None
                clinvar_negative = None
                lovd_positive = None
                lovd_negative = None
                for clinvar_sig in output_row["Clinical_Significance_ClinVar"].split(','):
                    if clinvar_sig in discordance_list_positive:
                        clinvar_positive = True
                    elif clinvar_sig in discordance_list_negative:
                        clinvar_negative = True
                for lovd_sig in output_row["Variant_effect_LOVD"].split(','):
                    if lovd_sig in discordance_list_positive:
                        lovd_positive = True
                    elif lovd_sig in discordance_list_negative:
                        lovd_negative = True
                if clinvar_positive is True and clinvar_negative is True:
                    output_row["Discordance_LOVD_And_ClinVar"] = True
                elif clinvar_positive is True and lovd_negative is True:
                    output_row["Discordance_LOVD_And_ClinVar"] = True
                elif lovd_positive is True and clinvar_negative is True:
                    output_row["Discordance_LOVD_And_ClinVar"] = True
                elif lovd_positive is True and lovd_negative is True:
                    output_row["Discordance_LOVD_And_ClinVar"] = True
                else:
                    output_row["Discordance_LOVD_And_ClinVar"] = False

            print output_row
            return output_row
        else:
            return False


    def handle(self, *args, **options):
        meta = {
            'file': '/tmp/discordance.csv',
            'class': CurrentVariant,
            'variant_fields': ('Genomic_Coordinate_hg38', 'HGVS_cDNA',
                       'Allele_frequency_ExAC', 'Allele_frequency_1000_Genomes'),
            'fields_of_interest': ('Genomic_Coordinate_hg38', 'HGVS_cDNA', 'Submitter_ClinVar',
                       'Clinical_Significance_ClinVar', 'Actionable_ClinVar',
                       'Discordance_ClinVar', 'Submitters_LOVD', 'Variant_effect_LOVD',
                       'Consistency_LOVD', 'Discordance_LOVD', 'Consistency_LOVD_And_ClinVar',
                       'Discordance_LOVD_And_ClinVar', 'Allele_frequency_ExAC',
                       'Allele_frequency_1000_Genomes')
        }

        f = open(meta['file'], 'w+')
        writer = csv.writer(f, encoding='utf-8')
        writer.writerow( meta['fields_of_interest'] )
        for obj in meta['class'].objects.all():
            obj_data = self._collect_data_for_variant(obj, meta)
            if not obj_data:
                continue
            else:
                writer.writerow(obj_data)
        f.close()
        print 'Data written to %s' % meta['file']

