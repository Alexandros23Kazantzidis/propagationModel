var script = document.createElement('script');
script.src = 'http://code.jquery.com/jquery-1.11.0.min.js';
script.type = 'text/javascript';
document.getElementsByTagName('head')[0].appendChild(script);

var body = document.getElementsByTagName("BODY")[0];

console.log(body);

body.onload = function () {

    var elem = document.getElementById("earth-model");
    console.log(elem.value);
};