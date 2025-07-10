import re
import plotly.graph_objects as graph_obj
from jinja2 import Template

def create_completeness_gauge_html(complete_percentages):
    categories = ['Essential Term Completeness', 'Pathway Coherence', 'Process Coherence', 'Protein Complex Coherence']
    completeness_data = dict(zip(categories, complete_percentages))

    html_template = Template("""
        <div class="gauge-table-container">
            <table class="gauge-table">
                <thead>
                    <tr>
                        <th>Category</th>
                        <th>Score</th>
                    </tr>
                </thead>
                <tbody>
                    {% for category, value in completeness_data.items() %}
                    <tr>
                        <td>{{ category }}</td>
                        <td>
                            <div id="gauge-{{ loop.index }}" class="gauge-container" style="width: 120px; height: 120px; margin: 0 auto;"></div>
                        </td>
                    </tr>
                    {% endfor %}
                </tbody>
            </table>
        </div>

        <script>
            document.addEventListener("DOMContentLoaded", function() {
                {% for category, value in completeness_data.items() %}
                Plotly.newPlot('gauge-{{ loop.index }}', [{
                    type: 'indicator',
                    mode: 'gauge+number',
                    value: {{ value }},
                    gauge: {
                        axis: { range: [0, 100] },
                        bar: { color: {{ "'#4caf50'" if value > 50 else "'#f44336'" }} },
                        steps: [
                            { range: [0, 50], color: '#f4f4f4' },
                            { range: [50, 100], color: '#e0e0e0' }
                        ]
                    }
                }], {
                    width: 120, height: 120, margin: { t: 0, b: 0, l: 0, r: 0 }
                });
                {% endfor %}
            });
        </script>
    """)

    return html_template.render(completeness_data=completeness_data)
