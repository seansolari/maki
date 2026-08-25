
import typer

from .extern_cli import app as real_data_app
from .synthetic_cli import app as synthetic_app

bmark_app = typer.Typer(help="Benchmarking module")


bmark_app.add_typer(real_data_app)
bmark_app.add_typer(synthetic_app, name="synthetic", help="Run synthetic benchmark")


@bmark_app.command("workflows", help="List available workflows")
def list_workflows():
    from .workflow import BenchmarkWorkflow
    
    for workflow in BenchmarkWorkflow.registry.values():
        typer.echo(f"{workflow.name}\t[{";".join(workflow.phase_names)}]\tsupports: {";".join(workflow.compatible_data_types)}")


if __name__ == "__main__":
    bmark_app()