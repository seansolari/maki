
import typer

from .extern_cli import app as real_data_app
from .synthetic_cli import app as synthetic_app

app = typer.Typer(help="Benchmarking module")


app.add_typer(real_data_app)
app.add_typer(synthetic_app, name="synthetic", help="Run synthetic benchmark")


@app.command("workflows", help="List available workflows")
def list_workflows():
    from .workflow import BenchmarkWorkflow
    
    for workflow in BenchmarkWorkflow.registry.values():
        typer.echo(f"{workflow.name}\t[{";".join(workflow.phase_names)}]\tsupports: {";".join(workflow.compatible_data_types)}")


if __name__ == "__main__":
    app()