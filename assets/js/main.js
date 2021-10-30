function render_plot(data_url, element_id) {
    return fetch(data_url)
        .then((resp) => resp.json())
        .then(function(fig) {
            canvas = document.getElementById(element_id);
            Plotly.newPlot(canvas, fig.data, fig.layout);
        })
}
