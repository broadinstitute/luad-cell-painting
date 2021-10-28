function render_plot(name) {
    return fetch(name + '.json')
        .then((resp) => resp.json())
        .then(function(fig) {
            canvas = document.getElementById(name);
            Plotly.newPlot(canvas, fig.data, fig.layout);
        })
}
