
import typer

from .data.extern.cli import app as real_data_app
from .data.synthetic.cli import app as synthetic_app
from .workflow import all_workflows

app = typer.Typer(help="Benchmarking module")


app.add_typer(real_data_app, name="dataset", help="Run benchmark on real dataset")
app.add_typer(synthetic_app, name="synthetic", help="Run synthetic benchmark")


@app.command("workflows", help="List available workflows")
def list_workflows():
    for workflow in all_workflows.values():
        typer.echo(f"{workflow.name}\t[{";".join(phase.name for phase in workflow.phases())}]\tsupports: {";".join(workflow.compatible_data_types)}")


if __name__ == "__main__":
    app()