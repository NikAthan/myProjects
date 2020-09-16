function click1(){
	var dots = document.getElementById("dots1");
	var moreText = document.getElementById("more1");
	var btnText = document.getElementById("myBtn");
		if (dots.style.display == "none") {
			dots.style.display = "inline";
			btnText.innerHTML = "Read more";
			btnText.style.borderWidth = "thin";
			moreText.style.display = "none";
		}
		else {
			dots.style.display = "none";
			btnText.innerHTML = "Read less";
			moreText.style.display = "inline";
		}
}

function click2(){
	var dots2 = document.getElementById("dots2");
	var moreText2 = document.getElementById("more2");
	var btn2Text = document.getElementById("Btn2");
		if (dots2.style.display == "none") {
			dots2.style.display = "inline";
			btn2Text.innerHTML = "Read more";
			moreText2.style.display = "none";
		}
		else {
			dots2.style.display = "none";
			btn2Text.innerHTML = "Read less";
			moreText2.style.display = "inline";
		}
}


$('.navbar-nav>li>a').on('click', function(){
    $('.navbar-collapse').collapse('hide');
});
