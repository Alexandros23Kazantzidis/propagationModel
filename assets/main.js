var script = document.createElement('script');
script.src = 'http://code.jquery.com/jquery-1.11.0.min.js';
script.type = 'text/javascript';
document.getElementsByTagName('head')[0].appendChild(script);

var body = document.getElementsByTagName("BODY")[0];

console.log(body);

window.onload = function () {

    console.log("break me");

    $('span#react-select-2--value-item').on('DOMSubtreeModified',function(){
        if ($('span[role="option"]')[0].textContent === "Keplerian") {
            $("#geodynamic-parameters").css("display", "none");
            $("#calculate-orbit").removeClass("disabled");
            $("#import-error").css("display","none");
        } else if ($('span[role="option"]')[0].textContent === "Geodynamic Model") {
            $("#geodynamic-parameters").css("display", "block");
            $("#calculate-orbit").addClass("disabled");
            $("#import-error").css("display","block");
        };
    });

    $('span#react-select-4--value-item').on('DOMSubtreeModified',function(){
        if (($('span[role="option"]')[2].textContent === "Harris-Priester")) {
            $("#goce").css("display", "none");
            $(".date-time").css("display", "block");
        } else if (($('span[role="option"]')[2].textContent === "Doorbnos Thermospheric Density Model (Only for GOCE)")) {
            $("#goce").css("display", "block");
            $(".date-time").css("display", "none");

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

    $('#calculate-geo-button').click(function() {
        $("#calculate-orbit").removeClass("disabled");
        $("#import-error").css("display","none");
    });

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

