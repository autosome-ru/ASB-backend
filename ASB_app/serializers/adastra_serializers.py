from flask_restplus import fields
from ASB_app.constants import chromosomes, bads, nucleotides


class ReleaseSerializers:
    def __init__(self, release):
        api = release.api
        self.release = release

        self.genome_polymorphism_location_model = api.model('Genome Mutation Position and id', {
            'chromosome': fields.String(enumerate=chromosomes),
            'position': fields.Integer,
            'ref': fields.String(enumerate=nucleotides),
            'alt': fields.String(enumerate=nucleotides),
            'rs_id': fields.Integer,
            'context': fields.String,
        })

        self.aggregated_snp_model = api.model('Aggregated SNP (no genome info) ', {
            'log_p_value_ref': fields.Float,
            'log_p_value_alt': fields.Float,
            'es_ref': fields.Float,
            'es_alt': fields.Float,
            'is_asb': fields.Boolean,
            'peak_calls': fields.Integer(min=0),
            'peak_callers': fields.Integer(min=0),
        })

        if int(self.release.version) >= 3:
            if int(self.release.version) >= 4:
                self.transcription_factor_model = api.model('Transcription factor', {
                    'tf_id': fields.Integer(readonly=True),
                    'name': fields.String,
                    'uniprot_ac': fields.String,
                    'gene_name': fields.String,
                    'aggregated_snps_count': fields.Integer(readonly=True),
                    'aggregated_snps_count005': fields.Integer(readonly=True),
                    'aggregated_snps_count010': fields.Integer(readonly=True),
                    'experiments_count': fields.Integer(readonly=True),
                })
            else:
                self.transcription_factor_model = api.model('Transcription factor', {
                    'tf_id': fields.Integer(readonly=True),
                    'name': fields.String,
                    'uniprot_ac': fields.String,
                    'aggregated_snps_count': fields.Integer(readonly=True),
                    'aggregated_snps_count005': fields.Integer(readonly=True),
                    'aggregated_snps_count010': fields.Integer(readonly=True),
                    'experiments_count': fields.Integer(readonly=True),
                })

            self.cell_line_model = api.model('Cell line', {
                'cl_id': fields.Integer(readonly=True),
                'name': fields.String,
                'aggregated_snps_count': fields.Integer(readonly=True),
                'aggregated_snps_count005': fields.Integer(readonly=True),
                'aggregated_snps_count010': fields.Integer(readonly=True),
                'experiments_count': fields.Integer(readonly=True),
            })
        else:
            self.transcription_factor_model = api.model('Transcription factor', {
                'tf_id': fields.Integer(readonly=True),
                'name': fields.String,
                'uniprot_ac': fields.String,
                'aggregated_snps_count': fields.Integer(readonly=True),
                'experiments_count': fields.Integer(readonly=True),
            })

            self.cell_line_model = api.model('Cell line', {
                'cl_id': fields.Integer(readonly=True),
                'name': fields.String,
                'aggregated_snps_count': fields.Integer(readonly=True),
                'experiments_count': fields.Integer(readonly=True),
            })

        if int(self.release.version) >= 3:
            self.gene_model = api.model('Gene', {
                'gene_id': fields.String(readonly=True),
                'gene_name': fields.String,
                'chromosome': fields.String(enumerate=chromosomes),
                'start_pos': fields.Integer,
                'end_pos': fields.Integer,
                'locus_start': fields.Integer,
                'locus_end': fields.Integer,
                'snps_count': fields.Integer,
                'snps_count010': fields.Integer,
                'eqtl_snps_count': fields.Integer,
                'eqtl_snps_count010': fields.Integer,
            })
        else:
            self.gene_model = api.model('Gene', {
                'gene_id': fields.String(readonly=True),
                'gene_name': fields.String,
                'chromosome': fields.String(enumerate=chromosomes),
                'start_pos': fields.Integer,
                'end_pos': fields.Integer,
                'locus_start': fields.Integer,
                'locus_end': fields.Integer,
                'snps_count': fields.Integer,
            })

        self.tf_snp_model = api.inherit('Transcription Factor SNP (no genome info)', self.aggregated_snp_model, {
            'tf_snp_id': fields.Integer(readonly=True),
            'motif_log_p_ref': fields.Float,
            'motif_log_p_alt': fields.Float,
            'motif_log_2_fc': fields.Float,
            'motif_orientation': fields.Boolean,
            'motif_position': fields.Integer,
            'motif_concordance': fields.String,
            'transcription_factor': fields.Nested(self.transcription_factor_model),
        })

        self.cl_snp_model = api.inherit('Cell line SNP (no genome info)', self.aggregated_snp_model, {
            'cl_snp_id': fields.Integer(readonly=True),
            'cell_line': fields.Nested(self.cell_line_model),
        })

        self.rs_snp_model = api.inherit('rs-SNP info for search', self.genome_polymorphism_location_model, {
            'dnase_aggregated_snps': fields.List(fields.Nested(self.cl_snp_model)),
            'atac_aggregated_snps': fields.List(fields.Nested(self.cl_snp_model)),
            'faire_aggregated_snps': fields.List(fields.Nested(self.cl_snp_model)),
            'has_concordance': fields.Boolean,
        })

        self.search_results_model = api.model('SNP search results model', {
            'total': fields.Integer,
            'results': fields.List(fields.Nested(self.rs_snp_model))
        })

        self.gene_search_results_model = api.model('SNP search results model, gene model included', {
            'total': fields.Integer,
            'gene': fields.Nested(self.gene_model),
            'results': fields.List(fields.Nested(self.rs_snp_model))
        })

        self.transcription_factor_browse_model = api.model('TF search results model', {
            'total': fields.Integer,
            'results': fields.List(fields.Nested(self.transcription_factor_model))
        })

        self.cell_line_browse_model = api.model('CellType search results model', {
            'total': fields.Integer,
            'results': fields.List(fields.Nested(self.cell_line_model))
        })

        self.exp_model_short = api.model('Experiment (short info)', {
            'exp_id': fields.String,
            'align': fields.String,
            'tf_name': fields.String(attribute='transcription_factor.name'),
            'cl_name': fields.String(attribute='cell_line.name'),
        })

        self.exp_snp_model = api.model('Experiment SNP', {
            'exp_snp_id': fields.Integer(readonly=True),
            'ref_readcount': fields.Integer,
            'alt_readcount': fields.Integer,
            'p_value_ref': fields.Float,
            'p_value_alt': fields.Float,
            'bad': fields.String(enumerate=bads),
            'experiment': fields.Nested(self.exp_model_short)
        })

        self.aggregated_snp_model_full = api.inherit('Aggregated SNP (with exp snps)', self.aggregated_snp_model, {
            'mean_bad': fields.Float,
            'exp_snps': fields.List(fields.Nested(self.exp_snp_model))
        })

        self.tf_snp_model_full = api.inherit('Transcription Factor SNP (with exp snps)', self.aggregated_snp_model_full, {
            'tf_snp_id': fields.Integer(readonly=True),
            'motif_log_p_ref': fields.Float,
            'motif_log_p_alt': fields.Float,
            'motif_log_2_fc': fields.Float,
            'motif_orientation': fields.Boolean,
            'motif_position': fields.Integer,
            'motif_concordance': fields.String,
            'transcription_factor': fields.Nested(self.transcription_factor_model),
        })

        self.cl_snp_model_full = api.inherit('Cell line SNP (with exp snps)', self.aggregated_snp_model_full, {
            'cl_snp_id': fields.Integer(readonly=True),
            'cell_line': fields.Nested(self.cell_line_model),
        })

        self.phenotype_model = api.model('Phenotype', {
            'db_name': fields.String,
            'phenotype_name': fields.String,
        })

        self.rs_snp_model_full = api.inherit('Complete rs-SNP info (with exp snps)', self.genome_polymorphism_location_model, {
            'atac_aggregated_snps': fields.List(fields.Nested(self.cl_snp_model_full)),
            'dnase_aggregated_snps': fields.List(fields.Nested(self.cl_snp_model_full)),
            'faire_aggregated_snps': fields.List(fields.Nested(self.cl_snp_model_full)),
            'phenotypes': fields.List(fields.Nested(self.phenotype_model))
        })

        if int(self.release.version) >= 3:
            self.frontpage_statistics_model = api.inherit('Front page statistics', {
                'transcription_factors_count': fields.Integer,
                'cell_types_count': fields.Integer,
                'snps_count': fields.Integer,
                'asbs_count': fields.Integer,
                'snps_count005': fields.Integer,
                'asbs_count005': fields.Integer,
                'snps_count010': fields.Integer,
                'asbs_count010': fields.Integer,
            })
        else:
            self.frontpage_statistics_model = api.inherit('Front page statistics', {
                'transcription_factors_count': fields.Integer,
                'cell_types_count': fields.Integer,
                'snps_count': fields.Integer,
            })
