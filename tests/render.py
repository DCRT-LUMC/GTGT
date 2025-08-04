import argparse
import json
import os
from typing import Any, Optional

from jinja2 import Environment, FileSystemLoader

payload = Optional[dict[str, Any]]


def render(
    variant: Optional[str], results: payload = None, links: payload = None
) -> str:
    template_file = "src/gtgt/templates/index.html.j2"
    template_folder = os.path.dirname(template_file)
    environment = Environment(loader=FileSystemLoader(template_folder))
    template = environment.get_template(os.path.basename(template_file))

    if not any((variant, results, links)):
        return template.render()
    elif variant and not any((results, links)):
        return template.render(variant=variant)
    elif variant and results:
        return template.render(variant=variant, results=results)
    raise RuntimeError()

    return template.render(varant=variant, results=results, links=links)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("results")
    parser.add_argument("--error", "-e", help="Error payload")

    args = parser.parse_args()

    with open(args.results) as fin:
        results = json.load(fin)

    print(render("10A>T", results))
