from datetime import datetime
import os


try:
    import pdf_reports.tools as pdf_tools
    from pdf_reports import (
        add_css_class,
        dataframe_to_html,
        pug_to_html,
        style_table_rows,
        write_report,
    )
    import pandas
    import matplotlib

    REPORT_PKGS_AVAILABLE = True
except ImportError:
    REPORT_PKGS_AVAILABLE = False


from .version import __version__

THIS_PATH = os.path.dirname(os.path.realpath(__file__))
ASSETS_PATH = os.path.join(THIS_PATH, "report_assets")
REPORT_TEMPLATE = os.path.join(ASSETS_PATH, "assessment_report.pug")
STYLESHEET = os.path.join(ASSETS_PATH, "report_style.css")


def end_pug_to_html(template, **context):
    now = datetime.now().strftime("%Y-%m-%d")
    defaults = {
        "sidebar_text": "Generated on %s by EGF's Plasmid assessor (version %s)"
        % (now, __version__),
        "plasma_logo_url": os.path.join(ASSETS_PATH, "imgs", "plasma_logo.png"),
    }
    for k in defaults:
        if k not in context:
            context[k] = defaults[k]
    return pug_to_html(template, **context)


def write_pdf_report(target, assessment):
    """Write an assessment report with a PDF summary.


    **Parameters**

    **target**
    > Path for PDF file.

    **assessment**
    > Assessment instance.
    """
    if not REPORT_PKGS_AVAILABLE:
        raise ImportError(
            "Install extra packages with `pip install plasmid_assessor[report]`"
        )

    assessment.figure_data = pdf_tools.figure_data(assessment.fig, fmt="svg")

    html = end_pug_to_html(REPORT_TEMPLATE, assessment=assessment,)
    write_report(html, target, extra_stylesheets=(STYLESHEET,))
