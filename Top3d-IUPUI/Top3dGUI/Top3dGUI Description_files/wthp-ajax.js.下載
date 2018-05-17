jQuery(document).ready(function($){
								
			 				
	$(".helpful-block-content a").click(function(event){
	
	  event.preventDefault();
	  
	  anc = $(this);
	  
	  main = anc.parent();
	 
	  href = anc.attr('href');
	 
	  main.find('.loader').remove();
	  
	  main.find('p').remove();
      
	  ajx_options = {
	  
		url: wthp.ajaxUrl+href,
		
		data: { action : 'wthp_post_helpful' },
		
		dataType: 'json',
		
		type : 'POST',
		
		success: function(resp){
			
		   elem = $('<p>').hide();
		   
				   
		   if( resp.success == true){
			 
			 elem.addClass('success');
		   
		   }else if( resp.error == true ){
			 
			 elem.addClass('error');
			 
		   }
		   
		   main.parent().html(elem);
		   
		   var c=elem.html(resp.message).fadeIn(500);
		   
		},
		
		beforeSend: function(){
		  
		  main.append('<div class="loader"></div>');
		  
		},
		
		complete: function(){
		
		 main.find('.loader').remove();	
		
		}
	  
	  };
	  
	   ajx_options.data.honeypot_spam = main.parent().find("input[name='honeypot']").val();
	    
	  if( anc.index() == 1 || anc.index() == 0 && wthp.positive_feedback == 1 ){
		 
		var mess_box = $('<textarea id="message" rows="3" cols="30"  style="margin-top:10px"></textarea>');
		
		var negative_btn = $('<input type="button" value="Submit"/>').css({ margin :"10px 3px 0", float: "right"});
         
        if(anc.index() == 1)
         {
		    main.parent().find('.wthp_title').text(wthp.hlp_title_nothank).parent().next().html(mess_box);
		 }
		else
		 {
			main.parent().find('.wthp_title').text(wthp.hlp_title_yesthank).parent().next().html(mess_box);
		 }

		mess_box.parent().append(negative_btn).append("<input type='hidden' value='"+href+"' />");
		
		mess_box.focus();
		
				
		
		negative_btn.click(function(){
		  
		  textarea_elem = $(this).parent().find("textarea");
	
		  border = mess_box.css('border');
	
		  mess_box.css('border', border);
		  
		 submit_negative_response(mess_box);  
		  
		})
		
		
		 return false;	
	  }
	  
	  
	  $.ajax(ajx_options);
	  
	
	})	
	
	function submit_negative_response(textarea){

    	message = $.trim(textarea.val()); 
			
		if( message ){
			   
			   main = textarea.parent();
			   
			   href = main.find("input[type='hidden']").val(); 
			   
			   ajx_options.url = wthp.ajaxUrl+href;
			   
			   ajx_options.data.message = message;	
			  
			   $.ajax(ajx_options);
        
			   
		}else{
			 
			  textarea.css( 'border', '1px solid red ');
				
		}
	
		
	}						
								
})

