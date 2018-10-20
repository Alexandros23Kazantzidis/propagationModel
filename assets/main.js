var script = document.createElement('script');
script.src = 'http://code.jquery.com/jquery-1.11.0.min.js';
script.type = 'text/javascript';
document.getElementsByTagName('head')[0].appendChild(script);

var body = document.getElementsByTagName("BODY")[0];

console.log(body);

window.onload = function () {

    $('span[role="option"]').on('DOMSubtreeModified',function(){
        if ($('span[role="option"]')[0].textContent === "Keplerian") {
            $("#geodynamic-parameters").css("display", "none");
        } else if ($('span[role="option"]')[0].textContent === "Geodynamic Model") {
            $("#geodynamic-parameters").css("display", "block");
        };
        if (($('span[role="option"]')[2].textContent === "Harris-Priester")) {
            $("#goce").css("display", "none");
        } else if (($('span[role="option"]')[2].textContent === "Doorbnos Thermospheric Density Model (Only for GOCE)")) {
            $("#goce").css("display", "block");
        }
    });

    jQuery('#other-forces').find('input').attr("class","other-forces-inputs");
    $('.other-forces-inputs')[0].className = "sun-moon-input";
    $('.other-forces-inputs')[0].className = "radiation-input";
    $('.other-forces-inputs')[0].className = "drag-input";

    $('.radiation-input').change(function () {
        if ($(this).is(":checked")) {
            $('#radiation-parameters').css('display', 'block');
        } else {
            $('#radiation-parameters').css('display', 'none');
        }
    });

    $('.drag-input').change(function () {
        if ($(this).is(":checked")) {
            $('#drag-parameters').css('display', 'block');
            $('#drag-model').css('display', 'block');
        } else {
            $('#drag-parameters').css('display', 'none');
             $('#drag-model').css('display', 'none');
        }
    });

    console.log($('#drag-model'));
//    $('.drag-input').change(function () {
//        if ($(this).is(":checked")) {
//            $('#drag-parameters').css('display', 'block');
//            $('#drag-model').css('display', 'block');
//        } else {
//            $('#drag-parameters').css('display', 'none');
//             $('#drag-model').css('display', 'none');
//        }
//    });


    jQuery(".Select-clear").hide();
};

